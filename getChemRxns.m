function [chems,S,rates,fast_species,fast_pairs,fast_affinity,S_cat, species_fast, stoic_fast, S_fast, S_cat_fast, r,p,r_fast,p_fast] = getChems(f)

    function [species,stoic,rates,r_species,r_stoic,p_species,p_stoic]=parseRxns(trans,fast)
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
                 rxn_trans = [chem_set  trans  chem_set '[^\,\n]*(\n|$)'];
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
            rates=cellfun(@(x) regexprep(x(end),'[\;\n\r]',''),rates,'UniformOutput',0);
            
            r_stoic=cf2(@(x)abs(x),lhs_stoic);
            p_stoic=cf2(@(x)abs(x),rhs_stoic);
            r_species=lhs_species;
            p_species=rhs_species;
            if strcmp(trans,rvsbl_rxn)
                
                rates=cellfun(@(x) strsplit(x{1},','),rates,'UniformOutput',0);
                rates=strtrim([rates{:}]);
                
                species=cellfun(@(x,y) [x,y],lhs_species,rhs_species,'UniformOutput',0);
                stoic=cellfun(@(x,y) [x,y],lhs_stoic,rhs_stoic,'UniformOutput',0);

                
                species=repelem(species,2);
                stoic=[stoic; cf2(@(x) -x,stoic)];
                stoic=stoic(:)';
%                 r_stoic=[r_stoic; cf2(@(x) -x,r_stoic)];
%                 r_stoic=r_stoic(:)';
%                 p_stoic=[p_stoic; cf2(@(x) -x,p_stoic)];
%                 p_stoic=p_stoic(:)';
                tmp=[r_stoic; p_stoic];
                p_stoic =[p_stoic; r_stoic];
                r_stoic=tmp(:)';
                p_stoic = p_stoic(:)';
                tmp = [r_species; p_species];
                p_species=[p_species; r_species];
                r_species = tmp(:)';
                p_species = p_species(:)';
%                 p_species = flipud(r_species);
            else
                nonempty = @(x) x(~cellfun(@isempty,x));
                rates=[rates{:}];
                rates=regexprep(rates,'[ ]+','');
                species=cellfun(@(l,r) [l r],lhs_species,rhs_species,'UniformOutput',false);
                stoic=cellfun(@(l,r) [l r],lhs_stoic,rhs_stoic,'UniformOutput',false);
            end
        else
            species={};
            stoic={};
            rates={};
            r_stoic={};
            p_stoic={};
            r_species={};
            p_species={};
        end
        
        
    end

str=fileread(f);

str=regexprep(str,"%[ \%\f\w\=\(\)\+\;\:\.\*\,\]\[\-\/\'\^\?]+",""); %remove comments
str=regexprep(str,"\'[^\'\n\r]+\'",""); %remove hardcoded strings with single quotes
str=regexprep(str,'\"[^\"\n\r]+\"',""); %remove hardcoded strings with double quotes
str=regexprep(str,'function[^\=]+\=[^\=]+\)',""); %remove function definition

name='[a-zA-Z_$][a-zA-Z_$0-9:\-]*';


str=regexprep(str,['[^\n]*D\(' name '\)[^\n]*\n'],""); %remove diffusion rate declarations
str=regexprep(str,['[^\n]*QSS(A?)\(' name '[^\n]*\n'],""); %remove QSS declarations

whitespace='[ \f\t\v]*';
nada={'0','{}'};
emptyset=strjoin(nada,'|');
empty_or_num='[ \f\t\v\.0-9]*';

chem_set = ['(' empty_or_num '(' name '|' emptyset ')' whitespace '[+])*?' '(' empty_or_num name '|' emptyset ')' whitespace ];

chem_set = ['((?:' empty_or_num '(' name '|' emptyset ')' whitespace '[+])*+(?:' empty_or_num name '|' emptyset '))' whitespace ];

forw_rxn = [  '[^<\na-zA-Z_$0-9]*-' whitespace '>' ];
backw_rxn = [ '<' whitespace '-[^>\na-zA-Z_$0-9]*'  ];
rvsbl_rxn = [ '<' whitespace '-' whitespace '>'];

trans = [whitespace '(' forw_rxn '|' backw_rxn '|' rvsbl_rxn '){1}+' whitespace];


rxn = [chem_set  trans  chem_set '[^\r\n]*'];

[rxns,chems]=regexp(str,rxn,"match","tokens");
chems=cellfun(@(x) x(~cellfun('isempty', x)),chems,'UniformOutput',0);
chems=regexp([chems{:}],['(' name ')'],'tokens');
chems=unique(cellfun(@(x) x(1),[chems{:}]),'stable');

species=cell(1,3);
stoic=cell(1,3);
p_stoic=cell(1,3);
r_stoic=cell(1,3);
p_species=cell(1,3);
r_species=cell(1,3);
rates=cell(1,3);

[species{1},stoic{1},rates{1},r_species{1},r_stoic{1},p_species{1},p_stoic{1}]=parseRxns(rvsbl_rxn);
[species{2},stoic{2},rates{2},r_species{2},r_stoic{2},p_species{2},p_stoic{2}]=parseRxns(forw_rxn);
[species{3},stoic{3},rates{3},r_species{3},r_stoic{3},p_species{3},p_stoic{3}]=parseRxns(backw_rxn);
%

species=[species{:}];
stoic=[stoic{:}];
rates=[rates{:}];
r_stoic=[r_stoic{:}];
p_stoic=[p_stoic{:}];
r_species=[r_species{:}];
p_species=[p_species{:}];

Nrx=size(species,2);
S=zeros(length(chems),Nrx);
r=zeros(length(chems),Nrx);
p=zeros(length(chems),Nrx);
S_cat=false(length(chems),Nrx);


for i=1:Nrx
    S(:,i)=cellfun(@(c) sum([stoic{i}{strcmp(c,species{i})}]),chems)';
    r(:,i)=cellfun(@(c) sum([r_stoic{i}{strcmp(c,r_species{i})}]),chems)';
    p(:,i)=cellfun(@(c) sum([p_stoic{i}{strcmp(c,p_species{i})}]),chems)';
    S_cat(:,i) = cellfun(@(c) ~isempty([stoic{i}{strcmp(c,species{i})}]),chems) == (0==S(:,i))';
end

if any(S(:)~=p(:)-r(:))
    error('Inconsistency in construction of stoichiometric matrix')
else
    S=p-r;
end



% S=[cell2mat(cellfun(@(x) ismember(chems,[x{:}]),lhs,'UniformOutput',0)');
% cell2mat(cellfun(@(x) ismember(chems,[x{:}]),rhs,'UniformOutput',0)')]

[species_fast,stoic_fast,affinity,r_species_fast,r_stoic_fast,p_species_fast, p_stoic_fast]=parseRxns(rvsbl_rxn,true);



Nrx_fast=size(species_fast,2);

S_fast=zeros(length(chems),Nrx_fast);
r_fast=zeros(length(chems),Nrx_fast);
p_fast=zeros(length(chems),Nrx_fast);
S_cat_fast=logical(S_fast);

for i=1:Nrx_fast
    S_fast(:,i)=cellfun(@(c) sum([stoic_fast{i}{strcmp(c,species_fast{i})}]),chems)';
    r_fast(:,i)=cellfun(@(c) sum([r_stoic_fast{i}{strcmp(c,r_species_fast{i})}]),chems)';
    p_fast(:,i)=cellfun(@(c) sum([p_stoic_fast{i}{strcmp(c,p_species_fast{i})}]),chems)';
    S_cat_fast(:,i) = cellfun(@(c) ~isempty([stoic_fast{i}{strcmp(c,species_fast{i})}]),chems) == (0==S_fast(:,i))';
end

if any(S(:)~=p(:)-r(:))
    error('Inconsistency in construction of stoichiometric matrix for fast reactions')
else
    S_fast=p_fast-r_fast;
end

fast_species=setdiff([species_fast{:}],[species{:}],'stable');

if any(abs(cell2mat([stoic_fast{:}]))~=1)

     error('only first-order instantaneous reactions are currently supported')

end

fast_pairs=cell(length(affinity),1);
fast_affinity=affinity;

% for i=1:2:length(species_fast)
%     try
%         fast_boy=any(reshape(cell2mat(cellfun(@(x) strcmp(fast_species,x),species_fast{i},'UniformOutput',false)),[length(fast_species),2]));
%     catch e
%         error(['Shape inconsistency in fast reaction processing.' newline 'Multiple reactions (unpaired) on the same line?'])
%     end
%     if nnz(fast_boy)>1
%         error('Intraconversion between fast species not currently supported')
%     end
%     j=ceil(i/2);
%     if fast_boy(1)
%         affinity{j}=['1/(' affinity{ceil(i/2)} ')'];
%     end
%     
%     fast_pairs{j}=species_fast{i}{~fast_boy};
%     fast_affinity{j}=affinity{ceil(i/2)};
%       
% 
% end


end

