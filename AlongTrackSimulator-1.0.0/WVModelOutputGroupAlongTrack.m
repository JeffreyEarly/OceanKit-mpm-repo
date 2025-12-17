classdef WVModelOutputGroupAlongTrack < WVModelOutputGroup
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        missionName
        ats

        repeatCycle
        tracks
        firstPassoverTime        
    end

    methods
        function self = WVModelOutputGroupAlongTrack(model,missionName,ats)
            arguments
                model WVModel
                missionName {mustBeText}
                ats AlongTrackSimulator
            end
            self@WVModelOutputGroup(model,name=missionName);
            self.missionName = missionName;
            self.ats = ats;
            self.repeatCycle = ats.repeatCycleForMissionWithName(missionName);

            if ~isinf(self.repeatCycle)
                wvt = self.model.wvt;
                alongtrack = ats.projectedPointsForRepeatMissionWithName(missionName,Lx=wvt.Lx,Ly=wvt.Ly,lat0=wvt.latitude,lon0=0);
                self.tracks = self.convertTrackVectorToPassoverCellArray(alongtrack);

                % only stop the model at the first t in the passover.
                self.firstPassoverTime = zeros(length(self.tracks),1);
                for iPassover = 1:length(self.tracks)
                    self.firstPassoverTime(iPassover) = self.tracks{iPassover}.t(1);
                end
            end
            self.addObservingSystem(WVAlongTrackObservingSystem(model,self));
        end

        function aString = description(self)
            if ~isinf(self.repeatCycle)
                aString = "Sampling pattern for the " + self.missionName + " mission with a repeat cycle of " + string(self.repeatCycle) + "s.";
            else
                aString = "Sampling pattern for the " + self.missionName + " mission.";
            end
        end

        function t = outputTimesForIntegrationPeriod(self,initialTime,finalTime)
            arguments (Input)
                self WVModelOutputGroup
                initialTime (1,1) double
                finalTime (1,1) double
            end
            arguments (Output)
                t (:,1) double
            end

            if ~isinf(self.repeatCycle)
                iCycle = floor(initialTime/self.repeatCycle);
                fCycle = floor(finalTime/self.repeatCycle);
                t = reshape(self.firstPassoverTime + (iCycle:fCycle)*self.repeatCycle,[],1);

                t(t<initialTime) = [];
                t(t>finalTime) = [];
                t = setdiff(t,self.timeOfLastIncrementWrittenToGroup);
            else
                wvt = self.model.wvt;
                % increment by 1s, typical along track sampling frequency
                t_full = ceil(initialTime):1:floor(finalTime);
                alongtrack = self.ats.projectedPointsForMissionWithName(self.missionName,Lx=wvt.Lx,Ly=wvt.Ly,lat0=wvt.latitude,lon0=0,time=t_full);
                self.tracks = self.convertTrackVectorToPassoverCellArray(alongtrack);

                self.firstPassoverTime = zeros(length(self.tracks),1);
                for iPassover = 1:length(self.tracks)
                    self.firstPassoverTime(iPassover) = self.tracks{iPassover}.t(1);
                end
                t = self.firstPassoverTime;
            end
        end

        function writeTimeStepToNetCDFFile(self,t)
            % Override the behavior of the superclass.
            % When we reach a time point where the model stops, we will
            % actually write all the time points for the passover. The
            % incrementsWrittenToGroup will accurately reflect the length
            % of the time dimension, but timeOfLastIncrementWrittenToGroup
            % will be the time at which the model stopped.
            if ( ~isempty(self.group) && t > self.timeOfLastIncrementWrittenToGroup )
                if ~isinf(self.repeatCycle)
                    iPassover = find( abs(self.firstPassoverTime - mod(t,self.repeatCycle)) < 1, 1, 'first');
                else
                    iPassover = find( abs(self.firstPassoverTime - t) < 1, 1, 'first');
                end
                firstOutputIndex = self.incrementsWrittenToGroup + 1;
                lastOutputIndex = firstOutputIndex + length(self.tracks{iPassover}.t) - 1;
                outputIndices = firstOutputIndex:lastOutputIndex;
                outputTimes = t + self.tracks{iPassover}.t - self.tracks{iPassover}.t(1);

                self.group.variableWithName('t').setValueAlongDimensionAtIndex(outputTimes,'t',outputIndices);

                for iObs = 1:length(self.observingSystems)
                    self.observingSystems(iObs).writeTimeStepToFile(self.group,outputIndices);
                end

                self.incrementsWrittenToGroup = outputIndices(end);
                self.timeOfLastIncrementWrittenToGroup = t(1);
            end
        end

    end

    methods (Static)

        function vars = classRequiredPropertyNames()
            vars = {'missionName','name','repeatCycle'};
        end

        function propertyAnnotations = classDefinedPropertyAnnotations()
            arguments (Output)
                propertyAnnotations CAPropertyAnnotation
            end
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);
            propertyAnnotations(end+1) = CAPropertyAnnotation('missionName','name the mission');
            propertyAnnotations(end+1) = CAPropertyAnnotation('name','name of output group');
            propertyAnnotations(end+1) = CANumericProperty('repeatCycle', {}, 's','orbital repeat cycle (Inf indicates that the orbit is non-repeat)');
            % propertyAnnotations(end+1) = CANumericProperty('firstPassoverTime', {}, 's','model time of first passover into the model domain');
        end

        function tracks = convertTrackVectorToPassoverCellArray(alongtrack)
            % alongtrack is a struct with fields (t,x,y)
            trackIndices = find(diff(alongtrack.t)>1);
            trackIndices(end+1) = length(alongtrack.t);
            startIndex = 1;
            tracks = cell(length(trackIndices),1);
            for i=1:length(trackIndices)
                endIndex = trackIndices(i);
                tracks{i}.x = alongtrack.x(startIndex:endIndex);
                tracks{i}.y = alongtrack.y(startIndex:endIndex);
                tracks{i}.t = alongtrack.t(startIndex:endIndex);
                startIndex = endIndex+1;
            end
        end
    end
end