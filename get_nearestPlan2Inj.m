function  [nearestPlan, nearestInj] = get_nearestPlan2Inj(InjPar)

% load ARA targets in MBA Inj Plan (n=381)
% stored as [381 x 1] arrays: Plan_injNo Plan_x Plan_y [nearestPlan, nearestInj] = get_nearestPlan2Inj(InjPar)Plan_z Plan_ara
load MBA_InjectionPlan_Targets.mat

%% get the nearest 2 Plan# for each injection
for i=1:numel(InjPar.brnID)
    d1x=Plan_x-InjPar.x(i);
    d1y=Plan_y-InjPar.y(i);
    d1z=Plan_z-InjPar.z(i);
    d1xyz=sqrt(d1x.^2+d1y.^2+d1z.^2);
    
    [u v]=sort(d1xyz);
    dd1(i,1:2)= u(1:2);
    ii1(i,1:2)= v(1:2);
    
    if isfield(InjPar,'DoubleInj')
        if InjPar.DoubleInj(i)
            d2x=Plan_x-InjPar.x2(i);
            d2y=Plan_y-InjPar.y2(i);
            d2z=Plan_z-InjPar.z2(i);
            d2xyz=sqrt(d2x.^2+d2y.^2+d2z.^2);
            
            [u v]=sort(d2xyz);
            dd2(i,1:2)= u(1:2);
            ii2(i,1:2)= v(1:2);
            
        else
            dd2(i,1:2)=NaN;
            ii2(i,1:2)=NaN;
        end;
    end
end;

% record the nearest 2 Plan# for each actual injection (and, when animal was double injected,every second injection too
nearestPlan.PlanNo1=Plan_injNo(ii1);
nearestPlan.d1=dd1;
if isfield(InjPar,'DoubleInj')
    iDouble=~isnan(ii2(:,1));
    nearestPlan.PlanNo2=NaN*nearestPlan.PlanNo1;
    nearestPlan.PlanNo2(iDouble,:)=Plan_injNo(ii2(iDouble,:));
    nearestPlan.d2=dd2;
end;

%% get the nearest 2 injections to each Plan#
for i=1:numel(Plan_injNo)
    
    xx=InjPar.x;
    yy=InjPar.y;
    zz=InjPar.z;
    ara_id=InjPar.ara_id;
    Ainj=InjPar.Ainj;
    brnID=InjPar.brnID;
    if isfield(InjPar,'DoubleInj')
        xx=[xx InjPar.x2];
        yy=[yy InjPar.y2];
        zz=[zz InjPar.z2];
        ara_id={ara_id{:} InjPar.ara_id2{:}};
        Ainj=[Ainj InjPar.Ainj2];
        brnID={InjPar.brnID{:} InjPar.brnID{:}};
    end;
    d2x=Plan_x(i)-xx;
    d2y=Plan_y(i)-yy;
    d2z=Plan_z(i)-zz;
    
    d2xyz=sqrt(d2x.^2+d2y.^2+d2z.^2);
    
    [u v]=sort(d2xyz);
    dd(i,1:2)= u(1:2);
    ii(i,1:2)= v(1:2);
end;
% something is wrong in identifying Ainj below ...
nearestInj.d=dd;
nearestInj.brnID=brnID(ii);
nearestInj.InjNo=Ainj(ii);
nearestInj.ara_id=ara_id(ii);
nearestInj.x=xx(ii);
nearestInj.y=yy(ii);
nearestInj.z=zz(ii);

