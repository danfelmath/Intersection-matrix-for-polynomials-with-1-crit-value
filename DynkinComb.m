function [Vout,IoutR]=DynkinComb(M,Ds)
%% USAGE:  
%%         [Vout,IoutR]=DynkinComb(vector_expression, int_expression); 
%%         M is partition of 9. Ds=0,1,2
%% PURPOSE: 
%%          Compute the possibles Dynkin diagrams associated to the 
%%          partition M of 9, according to Rules in Section 3.2 of 
%%          D. Lopez Garcia's thesis.
%% RETURN:  
%%          - IoutR: contains all Dynkin diagrams for the partition M
%%          IoutR(2,4,j)=1 if this Dynkin diagram satifies Rule 1.
%%          IoutR(2,5,j)=1 if this Dynkin diagram satifies Rule 2.
%%          - Vour: This array contains all the Dynkin diagrams which holds
%%          rules 1 and 2.
%% NOTE:    
%%           - If Ds=1, then the Dynkin diagram is the associated to
%%           Equation (3.9).
%%           - If Ds=2, then the Dynkin diagram is the associated to
%%           Equation (3.10).
%%           - If Ds=0, then the Dynkin diagram is the associated to
%%           the order [1 2 3; 4 5 6 ;7 8 9]' (see AVGZ)
%% EXAMPLE_: 
%%           [Vout,IoutR]=DynkinComb([4,2,2,1],0);
%%            Vout
if sum(M)~=9
    fprintf('M is not a partition of 9');
    Vout=zeros(3);
    IoutR=zeros(3);
else
N=length(M);
d=9;

%% 1) We construct the possible partition  of 1,2,...,d given by M
Lout=nchoosek(1:d,M(1));
for i=2:N
    Lout=PartitionD(Lout,M(i),d);
end
%% 2) We count the cases without redundancy in the values. For this, we choose
%% the minimum index of a set of identified cirtical values as the value, 
%% e.g: c_1=c_2=c_3=1, c_4=c_5=4.
Iout=zeros(1,d);
for i=1:size(Lout,1)
    Iout(i,Lout(i,1:M(1)))=min(Lout(i,1:M(1)));
   for j=2:N
       Vaux=Lout(i,sum(M(1:j-1))+1:sum(M(1:j))); % Element of Lout(i) in position M(j)
       Iout(i,Vaux)=min(Vaux); 
   end
end
%% We choose the Dynkni diagram of Eq 3.9, 3.10, or AVGZ basis
if Ds==1
    Il=[0 1 0;1 0 0;0 0 1]; Ir=[1 0 0;0 0 1;0 1 0];
elseif Ds==2
    Il=[1 0 0;0 0 1;0 1 0]; Ir=[1 0 0;0 0 1;0 1 0];
else 
    Il=eye(3); Ir=eye(3);
end

Iout=unique(Iout,'rows');
%% 3) We apply the rules in section 3.2 of D. Lopez Garcia's thesis. 
k=1; Ind1=[]; 
for i=1:size(Iout,1)
    Iaux=unique(Iout(i,1:d)); %% Values apeering in Iout
    Ind=0;
    %% Initializing the rules:
    a12=0; a45=0; a78=0; a13=0; a46=0; a79=0; a23=0; a56=0; a89=0;
    a14=0; a25=0; a36=0; a17=0; a28=0; a39=0; a47=0; a58=0; a69=0;
    a24=0; a1245=0; a27=0; a1278=0; a34=0; a1346=0; a35=0; a2356=0;
    a37=0; a1379=0; a38=0; a2389=0; a57=0; a4578=0; a67=0; a4679=0;
    a68=0; a5689=0;
    for j=1:N
        F=find(Iout(i,1:d)==Iaux(j));
        Ind(j,1:size(F,2))=F;
        a12=a12+prod(ismember([1,2],Ind(j,:)));
        a45=a45+prod(ismember([4,5],Ind(j,:))); 
        a78=a78+prod(ismember([7,8],Ind(j,:)));
        
        a13=a13+prod(ismember([1,3],Ind(j,:)));
        a46=a46+prod(ismember([4,6],Ind(j,:)));
        a79=a79+prod(ismember([7,9],Ind(j,:)));
        
        a23=a23+prod(ismember([2,3],Ind(j,:)));
        a56=a56+prod(ismember([5,6],Ind(j,:)));
        a89=a89+prod(ismember([8,9],Ind(j,:)));

        a14=a14+prod(ismember([1,4],Ind(j,:)));
        a25=a25+prod(ismember([2,5],Ind(j,:)));
        a36=a36+prod(ismember([3,6],Ind(j,:)));
        
        a17=a17+prod(ismember([1,7],Ind(j,:)));
        a28=a28+prod(ismember([2,8],Ind(j,:)));
        a39=a39+prod(ismember([3,9],Ind(j,:)));

        a47=a47+prod(ismember([4,7],Ind(j,:)));
        a58=a58+prod(ismember([5,8],Ind(j,:)));
        a69=a69+prod(ismember([6,9],Ind(j,:)));   
        
        a24=a24+prod(ismember([2,4],Ind(j,:)));
        a1245=a1245+prod(ismember([1,2,4,5],Ind(j,:)));
        a27=a27+prod(ismember([2,7],Ind(j,:)));
        a1278=a1278+prod(ismember([1,2,7,8],Ind(j,:)));
        
        a34=a34+prod(ismember([3,4],Ind(j,:)));
        a1346=a1346+prod(ismember([1,3,4,6],Ind(j,:)));
        a35=a35+prod(ismember([3,5],Ind(j,:)));
        a2356=a2356+prod(ismember([2,3,5,6],Ind(j,:)));
        a37=a37+prod(ismember([3,7],Ind(j,:)));
        a1379=a1379+prod(ismember([1,3,7,9],Ind(j,:)));
        a38=a38+prod(ismember([3,8],Ind(j,:)));
        a2389=a2389+prod(ismember([2,3,8,9],Ind(j,:)));
        
        a57=a57+prod(ismember([5,7],Ind(j,:)));
        a4578=a4578+prod(ismember([4,5,7,8],Ind(j,:)));
        
        a67=a67+prod(ismember([6,7],Ind(j,:)));
        a4679=a4679+prod(ismember([4,6,7,9],Ind(j,:)));
        a68=a68+prod(ismember([6,8],Ind(j,:)));
        a5689=a5689+prod(ismember([5,6,8,9],Ind(j,:)));
    end
    %% If any rule is zero then this Dynkin diagram is not valid
    %% Rule 1:
    R11=(a12==0 && a45==0 && a78==0) || (a12~=0 && a45~=0 && a78~=0);
    R12=(a13==0 && a46==0 && a79==0) || (a13~=0 && a46~=0 && a79~=0);
    R13=(a23==0 && a56==0 && a89==0) || (a23~=0 && a56~=0 && a89~=0);
    R14=(a14==0 && a25==0 && a36==0) || (a14~=0 && a25~=0 && a36~=0);
    R15=(a17==0 && a28==0 && a39==0) || (a17~=0 && a28~=0 && a39~=0);
    R16=(a47==0 && a58==0 && a69==0) || (a47~=0 && a58~=0 && a69~=0);
    
    %% Rule 2:
    R21=(a24==0)||(a1245~=0);
    R22=(a27==0)||(a1278~=0);
    R23=(a34==0)||(a1346~=0);
    R24=(a35==0)||(a2356~=0);
    R25=(a37==0)||(a1379~=0);
    R26=(a38==0)||(a2389~=0);
    R27=(a57==0)||(a4578~=0);
    R28=(a67==0)||(a4679~=0);
    R29=(a68==0)||(a5689~=0);
    %% Total Rule Rule 1 and Rule 2
    R1t=R11 && R12 && R13 && R14 && R15 && R16;
    R2t=R21 && R22 && R23 && R24 && R25 && R26 && R27 && R28 && R29;
    
    %%% For visualization
    IoutR(1:3,1:3,i)=Il*(reshape(Iout(i,1:d),[3,3]))*Ir;
    IoutR(2,4,i)=R1t;
    IoutR(2,5,i)=R2t;
    
    if R1t && R2t==1 %% Both rules are holded
        Ind1(k)=i;
        k=k+1;
    end
end
Ind1=nonzeros(Ind1);
if length(Ind1)>0
    for i=1:length(Ind1)
        Vout(1:3,1:3,i)=Il*(reshape(Iout(Ind1(i),1:d),[3,3]))*Ir;
    end
else
    fprintf('There is no Dynkin diagram holds the rules 1,2 for this partition M');
    Vout=zeros(3);
end
end





%%Auxiliary function
function Lout=PartitionD(Lin,N,d)
%% USAGE:  
%%         Lout=MonMatrix(matrix_expression,int_expression, int_expression)
%%         - Lin(j,:) is a vector with intputs between 1 and d.
%%         - N<=d
%% PURPOSE: 
%%          Compute the combinations for [1..d]\Lin(j,:) of N numbers and 
%%          add it to the initial row.
%% RETURN:  
%%          The rows in Lout are the vectors whose initial values are
%%          Lin(j,:) and the next N are the combination of [1..d]\Lin(j,:)
%%          for any j=1... #rows of Lin   
%% NOTE:    
%%           This function is used in DynkinComb.
%% EXAMPLE_: 
%%           Lin=[1,2];
%%           Lout=PartitionD(Lin,4,9)
S1=size(Lin,1);
S2=size(Lin,2);
In=1:d;
l=1;
Lout=zeros(1,S2+N);
for i=1:S1
    Dif=setdiff(In,Lin(i,:));
    Laux=nchoosek(Dif,N);
    for j=1:nchoosek(d-S2,N)
        Lout(l,:)=[Lin(i,:) Laux(j,:)];
        l=l+1;
    end
end
