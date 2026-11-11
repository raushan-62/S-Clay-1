/* S-CLAY1 browser material-point solver. Equations mirror s_clay1.py. */
'use strict';

const $ = (id) => document.getElementById(id);
const defaults = {lam:.20,kappa:.05,nu:.20,M:1.20,beta:.70,mu:50,p0:100,q0:0,pm0:250,v0:1.9,alpha0:0,etaK0:.35,protocol:'q3cu',evTotal:0,eqTotal:.30,steps:3000,substep:.0001,cycles:20,compareMCC:false};
let lastRun = null;
let chartData = null;
let plotProgress = 1;
let plotAnimationFrame = null;
let cancelPlotAnimation = null;

const getNumber = id => Number($(id).value);
const soilParams = () => ({lam:getNumber('lam'),kappa:getNumber('kappa'),nu:getNumber('nu'),M:getNumber('M'),beta:getNumber('beta'),mu:getNumber('mu'),pFloor:1e-6,alphaLimitRatio:.999,tol:1e-10,bisections:60});
const getState = () => ({v:getNumber('v0'),pm:getNumber('pm0'),alpha:getNumber('alpha0'),maxDrift:0,plasticMultiplier:0,substeps:0,epsPPl:0,epsQPl:0});
const fYield = (p,q,pm,a,M) => (q-a*p)**2-(M*M-a*a)*(pm-p)*p;
const gradYield = (p,q,pm,a,M) => { const x=q-a*p; return {fp:-2*a*x-(M*M-a*a)*(pm-2*p),fq:2*x,fpm:-(M*M-a*a)*p,fa:-2*p*x+2*a*(pm-p)*p}; };
const elastic = (p,q,v,de,prm) => {const K=v*p/prm.kappa,G=3*K*(1-2*prm.nu)/(2*(1+prm.nu));return {p:p+K*de[0],q:q+3*G*de[1],K,threeG:3*G};};

function validateParams(prm,sig,state){
  if(!(prm.lam>prm.kappa&&prm.kappa>0))throw Error('Soil parameters must satisfy λ > κ > 0.');
  if(!(prm.nu>-.99&&prm.nu<.49))throw Error('Poisson ratio ν′ must be between −0.99 and 0.49.');
  if(!(prm.M>0&&prm.beta>=0&&prm.mu>=0))throw Error('Require M > 0, β ≥ 0, and μ ≥ 0.');
  if(!(sig.p>prm.pFloor&&state.v>1&&state.pm>prm.pFloor&&Math.abs(state.alpha)<prm.alphaLimitRatio*prm.M))throw Error('Initial state must have p′ > 0, v > 1, p′m > 0 and |α| < M.');
  const f=fYield(sig.p,sig.q,state.pm,state.alpha,prm.M);
  if(f>prm.tol*Math.max(sig.p*sig.p,1))throw Error('Initial stress lies outside the initial yield surface. Increase p′m or adjust the initial stress / α.');
}
function valid(p,q,state,prm){
  if(![p,q,state.v,state.pm,state.alpha].every(Number.isFinite))throw Error('The integration produced a non-finite value. Reduce the strain increment or change parameters.');
  if(p<=prm.pFloor)throw Error(`Mean effective stress reached p′ = ${p.toPrecision(4)} kPa.`);
  if(state.v<=1||state.pm<=prm.pFloor)throw Error('Specific volume or yield size left its admissible range.');
  if(Math.abs(state.alpha)>=prm.alphaLimitRatio*prm.M)throw Error('Yield inclination approached |α| = M.');
}

function oneSubstep(point,de,prm){
  const {p,q,ev,eq,state}=point; const f0=fYield(p,q,state.pm,state.alpha,prm.M);
  if(f0>prm.tol*Math.max(p*p,1))throw Error('Current stress is outside its yield surface; the integration cannot continue.');
  const trial=elastic(p,q,state.v,de,prm); valid(trial.p,trial.q,state,prm);
  const ft=fYield(trial.p,trial.q,state.pm,state.alpha,prm.M);
  if(ft<=prm.tol*Math.max(trial.p*trial.p,1)){
    const next={p:trial.p,q:trial.q,ev:ev+de[0],eq:eq+de[1],state:{...state,v:state.v*Math.exp(-de[0])}};
    valid(next.p,next.q,next.state,prm); return next;
  }
  let lo=0,hi=1;
  for(let i=0;i<prm.bisections;i++){
    const m=(lo+hi)/2,pt=elastic(p,q,state.v,[de[0]*m,de[1]*m],prm);
    if(fYield(pt.p,pt.q,state.pm,state.alpha,prm.M)>0)hi=m;else lo=m;
  }
  const f=(lo+hi)/2,edge=elastic(p,q,state.v,[de[0]*f,de[1]*f],prm),rem=[de[0]*(1-f),de[1]*(1-f)];
  const g=gradYield(edge.p,edge.q,state.pm,state.alpha,prm.M), K=edge.K,threeG=edge.threeG,eta=edge.q/edge.p;
  const alphaRate=prm.mu*((.75*eta-state.alpha)*Math.max(g.fp,0)+prm.beta*(eta/3-state.alpha)*Math.abs(g.fq));
  const pmRate=state.pm*state.v/(prm.lam-prm.kappa)*g.fp;
  const H=-(g.fpm*pmRate+g.fa*alphaRate),den=g.fp*g.fp*K+g.fq*g.fq*threeG+H;
  const num=g.fp*K*rem[0]+g.fq*threeG*rem[1];
  if(!Number.isFinite(den)||den<=0)throw Error('Plastic consistency denominator is non-positive. Try smaller increments or review the parameter set.');
  let dl=num/den;
  if(dl < -prm.tol){const next={p:trial.p,q:trial.q,ev:ev+de[0],eq:eq+de[1],state:{...state,v:state.v*Math.exp(-de[0])}};valid(next.p,next.q,next.state,prm);return next;}
  dl=Math.max(0,dl);
  const p1=edge.p+K*(rem[0]-dl*g.fp),q1=edge.q+threeG*(rem[1]-dl*g.fq),a1=Math.max(-prm.alphaLimitRatio*prm.M,Math.min(prm.alphaLimitRatio*prm.M,state.alpha+dl*alphaRate));
  const pm1=state.pm+dl*pmRate,v1=state.v*Math.exp(-de[0]);
  const sUn={...state,v:v1,pm:pm1,alpha:a1,plasticMultiplier:state.plasticMultiplier+dl,epsPPl:(state.epsPPl||0)+dl*g.fp,epsQPl:(state.epsQPl||0)+dl*g.fq};valid(p1,q1,sUn,prm);
  const drift=Math.abs(fYield(p1,q1,pm1,a1,prm.M))/Math.max(p*p,1),denom=(prm.M*prm.M-a1*a1)*p1;
  if(denom<=0)throw Error('Yield ellipse became degenerate during drift correction.');
  const pmCorrected=p1+(q1-a1*p1)**2/denom;
  const next={p:p1,q:q1,ev:ev+de[0],eq:eq+de[1],state:{...sUn,pm:pmCorrected,maxDrift:Math.max(state.maxDrift,drift)}};
  valid(next.p,next.q,next.state,prm);return next;
}
function updateSubstepped(point,de,prm,maxStep){
  const n=Math.max(1,Math.ceil(Math.max(Math.abs(de[0]),Math.abs(de[1]))/maxStep));
  if(n>10000)throw Error(`This load step requires ${n} substeps; the maximum is 10,000. Increase the number of load steps.`);
  let out=point;const sub=[de[0]/n,de[1]/n];
  for(let i=0;i<n;i++)out=oneSubstep(out,sub,prm);
  out.state.substeps+=n;return out;
}
function getProtocolTotals(){
  const type=$('protocol').value;let ev=getNumber('evTotal'),eq=getNumber('eqTotal');
  if(type==='cu'||type==='q3cu'){ev=0;} else if(type==='extension'){ev=0;if(eq>0)eq=-eq;} else if(type==='cyclic'){ev=0;} else if(type==='iso'){eq=0;} else if(type==='k0'){eq=2*ev/3;}
  return {ev,eq};
}
function simulateProtocol(start,prm,tot,protocol,n,maxStep,cycles){
  const records=[{...start,state:{...start.state},fraction:0}];let point={...start,state:{...start.state}},criticalReached=false;
  if(protocol==='cyclic'){
    const amplitude=Math.abs(tot.eq),legs=3,total=cycles*legs*n;
    let index=0,current=0;
    for(let cycle=0;cycle<cycles;cycle++)for(const target of [amplitude,-amplitude,0]){
      const leg=(target-current)/n;
      for(let j=0;j<n;j++){
        point=updateSubstepped(point,[0,leg],prm,maxStep);index++;
        records.push({...point,state:{...point.state},fraction:index/total});
      }
      current=target;
    }
  }else {
    let stableCriticalSteps=0;
    for(let i=1;i<=n;i++){
      const previous=point.state;
      point=updateSubstepped(point,[tot.ev/n,tot.eq/n],prm,maxStep);
      records.push({...point,state:{...point.state},fraction:i/n});
      if(protocol==='q3cu'){
        const eta=point.q/point.p,dAlpha=Math.abs(point.state.alpha-previous.alpha),dPm=Math.abs(point.state.pm-previous.pm)/Math.max(Math.abs(previous.pm),1e-12);
        stableCriticalSteps=Math.abs(eta-prm.M)<1e-5&&dAlpha<1e-7&&dPm<1e-7?stableCriticalSteps+1:0;
        if(stableCriticalSteps>=10){criticalReached=true;break;}
      }
    }
  }
  return {records,criticalReached};
}
function createRun(){
  const prm=soilParams(),state=getState(),start={p:getNumber('p0'),q:getNumber('q0'),ev:0,eq:0,state};validateParams(prm,start,state);
  const n=Math.floor(getNumber('steps')),maxStep=getNumber('substep'),tot=getProtocolTotals(),protocol=$('protocol').value;
  if(!Number.isInteger(n)||n<1||n>20000)throw Error('Choose between 1 and 20,000 load steps.');
  if(!(maxStep>0&&Number.isFinite(maxStep)))throw Error('Internal step limit must be greater than zero.');
  const cycles=protocol==='cyclic'?Math.floor(getNumber('cycles')):0;
  if(protocol==='cyclic'&&(!Number.isInteger(cycles)||cycles<1||cycles>1000))throw Error('Choose between 1 and 1,000 cycles.');
  if(protocol==='cyclic'&&!(Math.abs(tot.eq)>0))throw Error('Cyclic deviatoric strain amplitude must be greater than zero.');
  if(protocol==='cyclic'&&cycles*3*n>100000)throw Error('Cyclic path has too many increments; reduce cycles or steps per leg.');
  const mainResult=simulateProtocol(start,prm,tot,protocol,n,maxStep,cycles),records=mainResult.records;
  let mccRecords=null,mccInitial=null,mccSizeAdjusted=false,mccError='';
  if($('compareMCC').checked){
    const pmNeeded=start.p+start.q*start.q/(prm.M*prm.M*start.p);
    const mccState={...state,pm:Math.max(state.pm,pmNeeded),alpha:0,maxDrift:0,plasticMultiplier:0,substeps:0,epsPPl:0,epsQPl:0};
    const mccStart={p:start.p,q:start.q,ev:0,eq:0,state:mccState},mccPrm={...prm,beta:0,mu:0};
    mccSizeAdjusted=mccState.pm>state.pm+1e-10;
    mccInitial={...mccStart,state:{...mccState}};
    try{validateParams(mccPrm,mccStart,mccState);mccRecords=simulateProtocol(mccStart,mccPrm,tot,protocol,n,maxStep,cycles).records;}
    catch(error){mccError=error.message;}
  }
  return {records,mccRecords,mccInitial,mccSizeAdjusted,mccError,prm,totals:tot,steps:records.length-1,protocol,criticalReached:mainResult.criticalReached,initial:{...start,state:{...state}}};
}

function updateInitialStatus(){
  const box=$('initialStatus');
  try{const prm=soilParams(),s=getState(),p=getNumber('p0'),q=getNumber('q0');validateParams(prm,{p,q},s);const f=fYield(p,q,s.pm,s.alpha,prm.M);box.className='inline-status';box.textContent=`Initial state admissible · f / p′₀² = ${(f/(p*p)).toExponential(2)}`;}
  catch(e){box.className='inline-status invalid';box.textContent=e.message;}
}
function syncProtocol(){
  const mode=$('protocol').value, note=$('protocolNote');
  if(mode==='cu'){ $('evTotal').value='0';$('eqTotal').value='.04';$('steps').value='100';note.textContent='Constant volume: Δεv = 0. Positive Δεq drives compression. Set the total deviatoric strain and load steps.';}
  else if(mode==='q3cu'){ $('evTotal').value='0';$('eqTotal').value='.30';$('steps').value='3000';$('substep').value='.0001';note.textContent='Assignment 6 Q3 CU path mapped to S-CLAY1: Δεv = 0, Δεq = 0.30 over up to 3,000 increments. Stops after 10 consecutive steps satisfy |q/p′ − M| < 1e−5, |Δα| < 1e−7, and |Δp′m/p′m| < 1e−7.';}
  else if(mode==='extension'){ $('evTotal').value='0';$('eqTotal').value='-.04';$('steps').value='100';note.textContent='Constant volume: Δεv = 0. Negative Δεq drives triaxial extension.';}
  else if(mode==='cyclic'){ $('evTotal').value='0';$('eqTotal').value='.0025';$('steps').value='20';note.textContent='Triangular undrained cycles: 0 → +amplitude → −amplitude → 0. The paper cautions that S-CLAY1 is unsuitable for cyclic loading; treat this as exploratory only.';}
  else if(mode==='iso'){ $('evTotal').value='.01';$('eqTotal').value='0';note.textContent='Isotropic strain path: deviatoric strain is zero. Positive volumetric strain is compression.';}
  else if(mode==='k0'){ $('evTotal').value='.02';$('eqTotal').value=(2*.02/3).toPrecision(5);note.textContent='Oedometer / K₀ approximation with dεr = 0, so dεq = 2 dεv / 3.';}
  else {note.textContent='Apply a fixed-ratio path between the entered total volumetric and deviatoric strains.';}
  $('evTotal').disabled=mode==='cu'||mode==='q3cu'||mode==='extension'||mode==='cyclic';$('eqTotal').disabled=mode==='iso'||mode==='k0';$('cyclesLabel').hidden=mode!=='cyclic';
  $('evLabel').style.opacity=$('evTotal').disabled?'.55':'1';$('eqLabel').style.opacity=$('eqTotal').disabled?'.55':'1';
}

const fmt=(n,d=3)=>Number(n).toLocaleString(undefined,{maximumFractionDigits:d,minimumFractionDigits:0});
function showResults(run){
  const end=run.records[run.records.length-1];$('metricP').innerHTML=`${fmt(end.p,2)}<small> kPa</small>`;$('metricQ').innerHTML=`${fmt(end.q,2)}<small> kPa</small>`;$('metricAlpha').textContent=fmt(end.state.alpha,4);$('metricDrift').textContent=end.state.maxDrift.toExponential(2);$('stepCount').textContent=`${end.state.substeps.toLocaleString()} internal substeps · ${run.steps} prescribed increments`;
  const note=$('comparisonNote');if(run.protocol==='q3cu'){const eta=end.q/end.p;note.textContent=`Question 3 data loaded into the S-CLAY1 formulation. Assumed S-CLAY1-only inputs: α₀ = 0, β = 0.70, μ = 50. Near-critical checks: |q/p′ − M| = ${Math.abs(eta-run.prm.M).toExponential(2)}; ${run.criticalReached?'criterion reached.':'maximum Δεq reached before criterion.'}`;note.hidden=false;}else if(run.mccRecords){const m=run.mccRecords.at(-1);note.textContent=`MCC final: p′ ${fmt(m.p,2)} kPa, q ${fmt(m.q,2)} kPa. The baseline uses α₀ = β = μ = 0${run.mccSizeAdjusted?`; its p′m₀ was increased to ${fmt(run.mccInitial.state.pm,2)} kPa so the same starting stress is admissible`:''}.`;}else if(run.mccError){note.textContent=`S-CLAY1 completed; MCC baseline stopped: ${run.mccError}`;}else note.textContent='MCC comparison overlay is off.';note.hidden=run.protocol!=='q3cu'&&!$('compareMCC').checked;
  $('csvBtn').disabled=false;$('pngBtn').disabled=true;chartData=run;$('porePressureCard').hidden=!['cu','q3cu','extension','cyclic'].includes(run.protocol);plotProgress=0;drawAll(run,0);
}

function niceRange(min,max){if(!Number.isFinite(min)||!Number.isFinite(max))return [0,1];if(min===max){const pad=Math.max(Math.abs(min)*.1,1e-4);return [min-pad,max+pad];}const pad=(max-min)*.08;return [min-pad,max+pad];}
function chart(canvasId,series,options={}){
  const canvas=$(canvasId),rect=canvas.getBoundingClientRect(),dpr=window.devicePixelRatio||1,w=Math.max(250,rect.width),h=Math.max(150,rect.height);
  canvas.width=Math.round(w*dpr);canvas.height=Math.round(h*dpr);const c=canvas.getContext('2d');c.scale(dpr,dpr);c.clearRect(0,0,w,h);
  const points=series.flatMap(s=>s.data||[]).filter(p=>Number.isFinite(p.x)&&Number.isFinite(p.y));
  if(!points.length){c.fillStyle='#9aa69f';c.font='11px Arial';c.fillText('Run a simulation to view the response',14,h/2);return;}
  let minX=Infinity,maxX=-Infinity,minY=Infinity,maxY=-Infinity;for(const p of points){minX=Math.min(minX,p.x);maxX=Math.max(maxX,p.x);minY=Math.min(minY,p.y);maxY=Math.max(maxY,p.y);}
  let [xmin,xmax]=options.xRange||niceRange(minX,maxX);
  let [ymin,ymax]=options.yRange||niceRange(minY,maxY);
  if(xmin===xmax)[xmin,xmax]=niceRange(xmin,xmax);if(ymin===ymax)[ymin,ymax]=niceRange(ymin,ymax);
  const L=47,R=12,T=10,B=29,PW=w-L-R,PH=h-T-B;
  const X=x=>L+(x-xmin)/(xmax-xmin)*PW,Y=y=>T+PH-(y-ymin)/(ymax-ymin)*PH;
  c.lineWidth=1;c.strokeStyle='#edf0ec';c.fillStyle='#8b9890';c.font='9px Consolas, monospace';c.textAlign='right';c.textBaseline='middle';
  const ticks=4;for(let i=0;i<=ticks;i++){const yy=T+PH*i/ticks,val=ymax-(ymax-ymin)*i/ticks;c.beginPath();c.moveTo(L,yy);c.lineTo(L+PW,yy);c.stroke();c.fillText(options.yFormat?options.yFormat(val):fmt(val,1),L-7,yy);}
  c.textAlign='center';c.textBaseline='top';for(let i=0;i<=ticks;i++){const xx=L+PW*i/ticks,val=xmin+(xmax-xmin)*i/ticks;c.beginPath();c.moveTo(xx,T);c.lineTo(xx,T+PH);c.stroke();c.fillText(options.xFormat?options.xFormat(val):fmt(val,1),xx,T+PH+7);}
  c.strokeStyle='#cdd6cf';c.beginPath();c.moveTo(L,T);c.lineTo(L,T+PH);c.lineTo(L+PW,T+PH);c.stroke();
  for(const s of series){let data=s.data||[];if(s.animate&&Number.isFinite(options.progress)){const last=Math.max(0,(data.length-1)*Math.max(0,Math.min(1,options.progress))),i=Math.floor(last);data=data.slice(0,i+1);if(i<s.data.length-1){const t=last-i,a=s.data[i],b=s.data[i+1];data.push({x:a.x+(b.x-a.x)*t,y:a.y+(b.y-a.y)*t});}}c.save();c.strokeStyle=s.color||'#267b64';c.lineWidth=s.width||1.8;c.setLineDash(s.dash||[]);c.globalAlpha=s.alpha??1;c.beginPath();let begun=false;for(const p of data){const x=X(p.x),y=Y(p.y);if(!Number.isFinite(x)||!Number.isFinite(y))continue;if(!begun){c.moveTo(x,y);begun=true;}else c.lineTo(x,y);}if(begun)c.stroke();if(s.marker&&data.length){const p=data[data.length-1];c.fillStyle=s.color||'#267b64';c.beginPath();c.arc(X(p.x),Y(p.y),3.2,0,Math.PI*2);c.fill();}c.restore();}
  c.fillStyle='#74847b';c.font='9px Arial';c.textAlign='center';c.textBaseline='bottom';c.fillText(options.xLabel||'',L+PW/2,h-1);c.save();c.translate(10,T+PH/2);c.rotate(-Math.PI/2);c.textBaseline='top';c.fillText(options.yLabel||'',0,0);c.restore();
}
function yieldCurve(pm,alpha,M,branches=1){const out=[];const end=Math.max(pm,1e-5);for(let i=0;i<=120;i++){const p=end*i/120,r=Math.max(0,(M*M-alpha*alpha)*(pm-p)*p),center=alpha*p;out.push({x:p,y:center+Math.sqrt(r)});}if(branches===2)for(let i=120;i>=0;i--){const p=end*i/120,r=Math.max(0,(M*M-alpha*alpha)*(pm-p)*p),center=alpha*p;out.push({x:p,y:center-Math.sqrt(r)});}return out;}
function drawAll(run,progress=1){
  const r=run.records,end=r[r.length-1],start=run.initial,M=run.prm.M;
  const m=run.mccRecords,overlay=m?[{data:m.map(a=>({x:a.eq*100,y:a.q})),color:'#80649c',width:1.6,dash:[6,4],animate:true}]:[];
  chart('chartQ',[{data:r.map(a=>({x:a.eq*100,y:a.q})),color:'#21745e',width:2,animate:true},...overlay],{progress,xLabel:'Deviatoric strain, εq [%]',yLabel:'Deviator stress, q [kPa]',xFormat:x=>fmt(x,2),yFormat:y=>fmt(y,1)});
  chart('chartPath',[{data:r.map(a=>({x:a.p,y:a.q})),color:'#21745e',width:2,animate:true},{data:yieldCurve(start.state.pm,start.state.alpha,M,2),color:'#96aaa0',width:1,dash:[4,4]},{data:yieldCurve(end.state.pm,end.state.alpha,M,2),color:'#d49a45',width:1.5},...(m?[{data:m.map(a=>({x:a.p,y:a.q})),color:'#80649c',width:1.6,dash:[6,4],animate:true},{data:yieldCurve(run.mccInitial.state.pm,0,M,2),color:'#80649c',width:1,dash:[2,3]}]:[])],{progress,xLabel:'Mean effective stress, p′ [kPa]',yLabel:'Deviator stress, q [kPa]',xFormat:x=>fmt(x,0),yFormat:y=>fmt(y,0)});
  chart('chartV',[{data:r.map(a=>({x:a.p,y:a.state.v-1})),color:'#21745e',width:2,animate:true},...(m?[{data:m.map(a=>({x:a.p,y:a.state.v-1})),color:'#80649c',width:1.6,dash:[6,4],animate:true}]:[])],{progress,xLabel:'Mean effective stress, p′ [kPa]',yLabel:'Void ratio, e',xFormat:x=>fmt(x,0),yFormat:y=>fmt(y,3)});
  chart('chartAlpha',[{data:r.map(a=>({x:a.eq*100,y:a.state.alpha})),color:'#21745e',width:2,animate:true},...(m?[{data:m.map(a=>({x:a.eq*100,y:a.state.alpha})),color:'#80649c',width:1.6,dash:[6,4],animate:true}]:[])],{progress,xLabel:'Deviatoric strain, εq [%]',yLabel:'Yield inclination, α',xFormat:x=>fmt(x,2),yFormat:y=>fmt(y,3)});
  chart('chartEPPl',[{data:r.map(a=>({x:a.eq*100,y:a.state.epsPPl||0})),color:'#21745e',width:2,animate:true}],{progress,xLabel:'Deviatoric strain, εq [%]',yLabel:'Plastic volumetric strain, εpᵖˡ',xFormat:x=>fmt(x,2),yFormat:y=>fmt(y,4)});
  chart('chartEQPl',[{data:r.map(a=>({x:a.eq*100,y:a.state.epsQPl||0})),color:'#21745e',width:2,animate:true}],{progress,xLabel:'Deviatoric strain, εq [%]',yLabel:'Plastic deviatoric strain, εqᵖˡ',xFormat:x=>fmt(x,2),yFormat:y=>fmt(y,4)});
  if(['cu','q3cu','extension','cyclic'].includes(run.protocol)){const p0=start.p,q0=start.q;chart('chartU',[{data:r.map(a=>({x:a.eq*100,y:(p0-a.p)+(a.q-q0)/3})),color:'#21745e',width:2,animate:true},...(m?[{data:m.map(a=>({x:a.eq*100,y:(p0-a.p)+(a.q-q0)/3})),color:'#80649c',width:1.6,dash:[6,4],animate:true}]:[])],{progress,xLabel:'Deviatoric strain, εq [%]',yLabel:'Estimated excess pore pressure, Δu [kPa]',xFormat:x=>fmt(x,2),yFormat:y=>fmt(y,1)});}
}

function animatePlots(run){
  if(plotAnimationFrame)cancelAnimationFrame(plotAnimationFrame);
  if(cancelPlotAnimation)cancelPlotAnimation();
  return new Promise(resolve=>{
    const duration=12000,start=performance.now();let lastFrame=0,finished=false;
    const finish=()=>{if(finished)return;finished=true;cancelPlotAnimation=null;plotAnimationFrame=null;resolve();};
    cancelPlotAnimation=finish;
    const frame=now=>{if(finished)return;if(now-lastFrame<40&&now-start<duration){plotAnimationFrame=requestAnimationFrame(frame);return;}lastFrame=now;plotProgress=Math.min(1,(now-start)/duration);drawAll(run,plotProgress);if(plotProgress<1){setStatus(`PLOTTING ${Math.round(plotProgress*100)}%`,'busy');plotAnimationFrame=requestAnimationFrame(frame);}else{$('pngBtn').disabled=false;setStatus('COMPLETE');finish();}};
    plotProgress=0;drawAll(run,0);plotAnimationFrame=requestAnimationFrame(frame);
  });
}

function downloadCsv(){if(!chartData)return;const head=['step','fraction','p_kPa','q_kPa','eps_v','eps_q','eps_p_plastic','eps_q_plastic','void_ratio_e','specific_volume_v','p_m_kPa','alpha','yield_f','normalized_pre_correction_drift','plastic_multiplier_cumulative','substeps','estimated_excess_u_kPa_constant_cell_pressure','mcc_p_kPa','mcc_q_kPa','mcc_v','mcc_p_m_kPa','mcc_yield_f'];const p0=chartData.initial.p,q0=chartData.initial.q;const rows=chartData.records.map((a,i)=>{const m=chartData.mccRecords?.[i];return [i,a.fraction,a.p,a.q,a.ev,a.eq,a.state.epsPPl||0,a.state.epsQPl||0,a.state.v-1,a.state.v,a.state.pm,a.state.alpha,fYield(a.p,a.q,a.state.pm,a.state.alpha,chartData.prm.M),a.state.maxDrift,a.state.plasticMultiplier,a.state.substeps,(p0-a.p)+(a.q-q0)/3,m?.p??'',m?.q??'',m?.state.v??'',m?.state.pm??'',m?fYield(m.p,m.q,m.state.pm,m.state.alpha,chartData.prm.M):''];});const csv=[head,...rows].map(row=>row.map(x=>typeof x==='number'?Number(x).toPrecision(12):x).join(',')).join('\r\n');downloadBlob(new Blob([csv],{type:'text/csv;charset=utf-8'}),'s-clay1-response.csv');}
function downloadBlob(blob,name){const a=document.createElement('a'),url=URL.createObjectURL(blob);a.href=url;a.download=name;document.body.append(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(url),1000);}
function downloadCharts(){const ids=['chartQ','chartPath','chartV','chartAlpha','chartEPPl','chartEQPl'];if(['cu','q3cu','extension','cyclic'].includes(chartData.protocol))ids.push('chartU');const rows=Math.ceil(ids.length/2),w=1400,h=120+rows*380,out=document.createElement('canvas');out.width=w;out.height=h;const c=out.getContext('2d');c.fillStyle='#f4f5f0';c.fillRect(0,0,w,h);c.fillStyle='#152a25';c.font='700 26px Arial';c.fillText('S-CLAY1 · material point response',48,48);c.fillStyle='#74817a';c.font='13px Consolas';c.fillText(`p′=${fmt(chartData.records.at(-1).p,2)} kPa · q=${fmt(chartData.records.at(-1).q,2)} kPa · α=${fmt(chartData.records.at(-1).state.alpha,4)}`,48,74);ids.forEach((id,i)=>{const src=$(id),x=40+(i%2)*680,y=100+Math.floor(i/2)*380;c.fillStyle='#fff';c.fillRect(x,y,640,350);c.drawImage(src,x+10,y+25,620,310);});out.toBlob(b=>downloadBlob(b,'s-clay1-charts.png'));}

function setStatus(label,kind='ready'){$('modelStatus').innerHTML=`<i></i> ${label}`;$('modelStatus').style.color=kind==='error'?'#a84e43':kind==='busy'?'#ae7b35':'#37765f';$('modelStatus').querySelector('i').style.background=kind==='error'?'#c46f62':kind==='busy'?'#d49a45':'#42a77c';}
async function runSimulation(){const btn=$('runBtn'),err=$('errorBox');btn.disabled=true;err.hidden=true;setStatus('CALCULATING','busy');await new Promise(resolve=>requestAnimationFrame(()=>setTimeout(resolve,30)));try{const run=createRun();lastRun=run;showResults(run);await animatePlots(run);}catch(e){err.textContent=e.message;err.hidden=false;setStatus('CHECK INPUTS','error');}finally{btn.disabled=false;}}
function reset(){if(cancelPlotAnimation)cancelPlotAnimation();if(plotAnimationFrame)cancelAnimationFrame(plotAnimationFrame);plotProgress=1;for(const [id,val] of Object.entries(defaults)){if($(id)){if(id==='compareMCC')$(id).checked=val;else $(id).value=val;}}syncProtocol();lastRun=null;chartData=null;$('porePressureCard').hidden=true;$('comparisonNote').hidden=true;$('errorBox').hidden=true;$('csvBtn').disabled=true;$('pngBtn').disabled=true;$('metricP').innerHTML='—<small> kPa</small>';$('metricQ').innerHTML='—<small> kPa</small>';$('metricAlpha').textContent='—';$('metricDrift').textContent='—';$('stepCount').textContent='No run yet';setStatus('READY');updateInitialStatus();drawEmpty();}
function drawEmpty(){for(const id of ['chartQ','chartPath','chartV','chartAlpha','chartEPPl','chartEQPl','chartU'])chart(id,[{data:[]}],{xLabel:'',yLabel:''});}

$('runBtn').addEventListener('click',runSimulation);$('resetBtn').addEventListener('click',reset);$('csvBtn').addEventListener('click',downloadCsv);$('pngBtn').addEventListener('click',downloadCharts);$('protocol').addEventListener('change',()=>{syncProtocol();updateInitialStatus();});$('alphaK0Btn').addEventListener('click',()=>{const M=getNumber('M'),eta=getNumber('etaK0');if(!(M>Math.abs(eta))){$('initialStatus').textContent='Need M > |ηK₀| to calculate αK₀.';$('initialStatus').className='inline-status invalid';return;}const a=eta-.75*(M*M-eta*eta);$('alpha0').value=Number(a.toFixed(5));updateInitialStatus();});
$('q3PresetBtn').addEventListener('click',()=>{const preset={lam:.20,kappa:.05,nu:.20,M:1.20,beta:.70,mu:50,p0:100,q0:0,pm0:250,v0:1.90,alpha0:0,protocol:'q3cu',compareMCC:false};for(const [id,value] of Object.entries(preset)){if(id==='compareMCC')$(id).checked=value;else $(id).value=value;}syncProtocol();updateInitialStatus();$('comparisonNote').hidden=true;});
$('yieldSizeBtn').addEventListener('click',()=>{const p=getNumber('p0'),q=getNumber('q0'),a=getNumber('alpha0'),M=getNumber('M'),den=(M*M-a*a)*p;if(!(p>0&&den>0)){ $('initialStatus').textContent='Need p′ > 0 and |α| < M to set the yield size.';$('initialStatus').className='inline-status invalid';return;}const pm=p+(q-a*p)**2/den;$('pm0').value=Number(pm.toFixed(5));updateInitialStatus();});
for(const el of document.querySelectorAll('input'))el.addEventListener('input',updateInitialStatus);
window.addEventListener('resize',()=>{if(chartData)drawAll(chartData,plotProgress);else drawEmpty();});
syncProtocol();updateInitialStatus();drawEmpty();

