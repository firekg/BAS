function PL = a_Combine_Data(varargin)

fields = fieldnames(varargin{1});
for i = 1:numel(fields)
    PL.(fields{i})=[];
    for j=1:length(varargin)
      PL.(fields{i})=[PL.(fields{i});varargin{j}.(fields{i})];
    end
end

PL.Trials = 0;
for j = 1:length(varargin)
    PL.Trials = PL.Trials+varargin{j}.Trials;
end