function [chems,S,rates,fast_species,fast_pairs,fast_affinity] = getChems(f)

    function [species,stoic,rates]=parseRxns(trans,fast)
        if nargin==1
            fast=false;
        end
        
        cf2=@(f,x)  cellfun( @(y) cellfun(@(z) f(z),y,'UniformOutput',0),x,'UniformOutput',0);
        
        if strcmp(trans,rvsbl_rxn)||strcmp(trans,forw_rxn)
            f_lhs = -1;
            f_rhs = 1;
        elseif strcmp(trans,backw_rxn)
            f_lhs = 1;
            f_rhs = -1;
        end
        
        if strcmp(trans,rvsbl_rxn)
            if ~fast
                rxn_trans = [chem_set  trans  chem_set '[^\r\n]*\,[^\r\n]*'];
            else
                 rxn_trans = [chem_set  trans  chem_set '[^\,\n]*\n'];
            end
        else
            rxn_trans = [chem_set  trans  chem_set '[^\r\n]*'];
        end
        [rxns_trans]=regexp(str,rxn_trans,"match");
        if ~isempty(rxns_trans)
            
            rxn_split=regexp(rxns_trans,trans,"split");
            
            lhs=cellfun(@(x) regexp(x(1),['(' empty_or_num ')(' name ')'],'tokens'),rxn_split);
            %lhs_split=cf2(@(x) regexp(x,empty,"split"),lhs)                          %cellfun( @(x) cellfun(@(y) regexp(y,empty,"split"),x,'UniformOutput',0),lhs,'UniformOutput',0);
            
            lhs_stoic=cf2(@(x) f_lhs*empty2one(str2double(x{1})),lhs);                   %cellfun(@(x) cellfun(@(y) empty2one(str2double(y{1})),x,'UniformOutput',0),lhs,'UniformOutput',0);
            lhs_species=cf2(@(x) x{2},lhs);                                         %cellfun(@(x) cellfun(@(y) y{2},x,'UniformOutput',0),lhs,'UniformOutput',0);
            
            rhs=cellfun(@(x) regexp(x(2),chem_set,"match",'once'),rxn_split);
            rhs=cellfun(@(x) regexp(x,['(' empty_or_num ')(' name ')'],'tokens'),rhs,'UniformOutput',0);
            
            rhs_stoic=cf2(@(x) f_rhs*empty2one(str2double(x{1})),rhs);                %cellfun(@(x) cellfun(@(y) empty2one(str2double(y{1})),x,'UniformOutput',0),rhs,'UniformOutput',0)
            rhs_species=cf2(@(x) x{2},rhs);                                       %cellfun(@(x) cellfun(@(y) y{2},x,'UniformOutput',0),rhs,'UniformOutput',0)
            
            rates=cellfun(@(x) regexp(x(2),chem_set,"split",'once'),rxn_split);
            
            if strcmp(trans,rvsbl_rxn)
                
                rates=cellfun(@(x) regexprep(x(end),'[\;\n\r]',''),rates,'UniformOutput',0);
                rates=cellfun(@(x) strsplit(x{1},','),rates,'UniformOutput',0);
                rates=strtrim([rates{:}]);
                
                species=cellfun(@(x,y) [x,y],lhs_species,rhs_species,'UniformOutput',0);
                stoic=cellfun(@(x,y) [x,y],lhs_stoic,rhs_stoic,'UniformOutput',0);
                
                species=repelem(species,2);
                stoic=[stoic; cf2(@(x) -x,stoic)];
                stoic=stoic(:)';
            else
                species=[lhs_species, rhs_species];
                stoic=[lhs_stoic rhs_stoic];
            end
        else
            species={};
            stoic={};
            rates={};
        end
        
        
    end

str=fileread(f);

str=regexprep(str,"%[ \%\f\w\=\(\)\+\;\:\.\*\,\]\[\-\/\'\^]+",""); %remove comments
str=regexprep(str,"\'[^\'\n\r]+\'",""); %remove hardcoded strings with single quotes
str=regexprep(str,'\"[^\"\n\r]+\"',""); %remove hardcoded strings with double quotes
str=regexprep(str,'function[^\=]+\=[^\=]+\)',""); %remove function definition

name='[a-zA-Z_$][a-zA-Z_$0-9\-]*';
empty='[ \f\t\v]*';
empty_or_num='[ \f\t\v\.0-9]*';

chem_set = ['(' empty_or_num '(' name ')' empty '[+])*' '(' empty_or_num name ')' empty ];


forw_rxn = [  '[^<]*-' empty '>' ];
backw_rxn = [ '<' empty '-[^>]*'  ];
rvsbl_rxn = [ '<' empty '-' empty '>'];

trans = [empty '(' forw_rxn '|' backw_rxn '|' rvsbl_rxn ')' empty];


rxn = [chem_set  trans  chem_set '[^\r\n]*'];

[rxns,chems]=regexp(str,rxn,"match","tokens");
chems=cellfun(@(x) x(~cellfun('isempty', x)),chems,'UniformOutput',0);
chems=regexp([chems{:}],['(' name ')'],'tokens');
chems=unique(cellfun(@(x) x(1),[chems{:}]),'stable');

species=cell(1,3);
stoic=cell(1,3);
rates=cell(1,3);

[species{1},stoic{1},rates{1}]=parseRxns(rvsbl_rxn);
[species{2},stoic{2},rates{2}]=parseRxns(forw_rxn);
[species{3},stoic{3},rates{3}]=parseRxns(backw_rxn);
%

species=[species{:}];
stoic=[stoic{:}];
rates=[rates{:}];

Nrx=size(species,2);
S=zeros(length(chems),Nrx);



for i=1:Nrx
    S(ismember(chems,[species{i}]),i)=[stoic{i}{:}]';
end

% S=[cell2mat(cellfun(@(x) ismember(chems,[x{:}]),lhs,'UniformOutput',0)');
% cell2mat(cellfun(@(x) ismember(chems,[x{:}]),rhs,'UniformOutput',0)')]

[species_fast,stoic_fast,affinity]=parseRxns(rvsbl_rxn,true);

fast_species=setdiff([species_fast{:}],[species{:}]);

if any(abs(cell2mat([stoic_fast{:}]))~=1)

     error('only unimolecular instantaneous reactions are currently supported')

end

fast_pairs=cell(length(species_fast),1);
fast_affinity=cell(length(species_fast),1);

for i=1:2:length(species_fast)
    fast_boy=any(reshape(cell2mat(cellfun(@(x) strcmp(fast_species,x),species_fast{i},'UniformOutput',false)),[length(fast_species),2]));
    if nnz(fast_boy)>1
        error('intraconversion between fast species not currently supported')
    end
    if fast_boy(1)
        affinity{ceil(i/2)}=['1/(' affinity{ceil(i/2)} ')'];
    end
    
    fast_pairs{i}=species_fast{i}{~fast_boy};
    fast_affinity{i}=affinity{ceil(i/2)};
      

end


end

