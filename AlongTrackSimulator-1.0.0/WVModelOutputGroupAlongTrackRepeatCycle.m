classdef WVModelOutputGroupAlongTrackRepeatCycle < WVModelOutputGroup
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        missionName
        tracks
        repeatCycle

        firstPassoverTime
    end

    methods
        function self = WVModelOutputGroupAlongTrackRepeatCycle(model,missionName,tracks,repeatCycle)
            arguments
                model WVModel
                missionName {mustBeText}
                tracks
                repeatCycle
            end
            self@WVModelOutputGroup(model,missionName);
            self.missionName = missionName;
            self.tracks = tracks;
            self.repeatCycle = repeatCycle;

            % We treat the entire passover as one instant in time, and thus
            % only stop the model at the first t in the passover.
            self.firstPassoverTime = zeros(length(self.tracks),1);
            for iPassover = 1:length(self.tracks)
                self.firstPassoverTime(iPassover) = self.tracks{iPassover}.t(1);
            end

            self.addObservingSystem(WVAlongTrackObservingSystem(model,missionName,tracks,repeatCycle));
        end

        function aString = description(self)
            aString = "Sampling pattern for the " + self.missionName + " mission with a repeat cycle of " + string(self.repeatCycle) + "s.";
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
            iCycle = floor(initialTime/self.repeatCycle);
            fCycle = floor(finalTime/self.repeatCycle);
            t = reshape(self.firstPassoverTime + (iCycle:fCycle)*self.repeatCycle,[],1);
            
            t(t<initialTime) = [];
            t(t>finalTime) = [];
            t = setdiff(t,self.timeOfLastIncrementWrittenToGroup);
        end

        function writeTimeStepToNetCDFFile(self,t)
            % Override the behavior of the superclass.
            % When we reach a time point where the model stops, we will
            % actually write all the time points for the passover. The
            % incrementsWrittenToGroup will accurately reflect the length
            % of the time dimension, but timeOfLastIncrementWrittenToGroup
            % will be the time at which the model stopped.
            if ( ~isempty(self.group) && t > self.timeOfLastIncrementWrittenToGroup )
                iPassover = find( abs(self.firstPassoverTime - mod(t,self.repeatCycle)) < 1, 1, 'first');
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
end