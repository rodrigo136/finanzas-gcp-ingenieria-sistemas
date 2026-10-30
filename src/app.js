(function(){
  const LA=window.LinearAlgebra;
  const units=["Desarrollo","Infraestructura","Seguridad","Soporte"];
  const scenarios={
    ideal:{name:"Escenario ideal",desc:"Situación financiera estable. La matriz está equilibrada y bien condicionada; el vector b representa metas de asignación habituales.",A:[[10,2,1,.5],[2,9,1.5,.5],[1,1.5,8,1],[.5,.5,1,7]],b:[120,110,95,80]},
    stress:{name:"Escenario bajo estrés financiero",desc:"Se mantienen las relaciones estructurales entre unidades y se incrementan las metas ponderadas para representar inflación, mayores costos de servicios tecnológicos y presión operativa.",A:[[10,2,1,.5],[2,9,1.5,.5],[1,1.5,8,1],[.5,.5,1,7]],b:[150,138,125,105]},
    ill:{name:"Escenario mal condicionado",desc:"Se modelan escalas muy diferentes de ponderación entre unidades y acoplamientos relativos. La matriz se construye como D·C·D, con D=diag(1,10,100,1000) y C simétrica definida positiva. Esto produce un fuerte mal escalamiento diagonal, para el cual Jacobi es el precondicionador pertinente.",A:[[4,6,20,100],[6,100,200,1000],[20,200,40000,10000],[100,1000,10000,4000000]],b:[2170,21680,955600,60228000]}
  };
  let current="ideal", last=null;
  const $=id=>document.getElementById(id);
  const fmt=(v,d=6)=>Number.isFinite(v)?Number(v).toLocaleString("es-BO",{maximumFractionDigits:d}):"∞";
  const sci=v=>Number.isFinite(v)?v.toExponential(3):"∞";

  function renderEditors(){
    const A=scenarios[current].A,b=scenarios[current].b;
    $("matrixEditor").innerHTML=A.map((row,i)=>row.map((v,j)=>`<label class="sr-only" for="a${i}${j}">a${i+1}${j+1}</label><input id="a${i}${j}" type="number" step="any" value="${v}" aria-label="Coeficiente fila ${i+1}, columna ${j+1}">`).join("")).join("");
    $("vectorEditor").innerHTML=b.map((v,i)=>`<label class="vector-cell"><span>b${i+1}</span><input id="b${i}" type="number" step="any" value="${v}" aria-label="Término independiente ${i+1}"></label>`).join("");
    $("scenarioDescription").innerHTML=`<b>${scenarios[current].name}.</b> ${scenarios[current].desc}`;
    updateDiagnostics();
  }
  function readData(){
    const A=Array.from({length:4},(_,i)=>Array.from({length:4},(_,j)=>Number($(`a${i}${j}`).value)));
    const b=Array.from({length:4},(_,i)=>Number($(`b${i}`).value));
    if(A.flat().concat(b).some(v=>!Number.isFinite(v)))throw new Error("Todos los coeficientes deben ser números finitos.");
    return {A,b};
  }
  function getOptions(){
    const tol=Number($("tolerance").value), maxIter=Number($("maxIterations").value);
    if(!(tol>0&&tol<1))throw new Error("La tolerancia debe estar entre 0 y 1.");
    if(!(Number.isInteger(maxIter)&&maxIter>0&&maxIter<=1000))throw new Error("Las iteraciones máximas deben ser un entero entre 1 y 1000.");
    return {tol,maxIter,x0:$("initialVector").value==="ones"?[1,1,1,1]:[0,0,0,0]};
  }
  function updateDiagnostics(){
    try{
      const {A}=readData(), sym=LA.isSymmetric(A);
      $("symmetryMetric").textContent=sym?"Sí":"No";
      if(sym){LA.cholesky(A);$("spdMetric").textContent="Sí";$("spdBadge").textContent="Matriz válida";$("spdBadge").className="status-pill good";$("conditionMetric").textContent=fmt(LA.conditionNumber(A),2);$("preConditionMetric").textContent=fmt(LA.preconditionedCondition(A),2);}
      else{$("spdMetric").textContent="No";$("spdBadge").textContent="No simétrica";$("spdBadge").className="status-pill bad";$("conditionMetric").textContent="—";$("preConditionMetric").textContent="—";}
    }catch(e){$("spdMetric").textContent="No";$("spdBadge").textContent="No SPD";$("spdBadge").className="status-pill bad";$("conditionMetric").textContent="—";$("preConditionMetric").textContent="—";$("symmetryMetric").textContent="—";}
  }
  function run(){
    try{
      const {A,b}=readData(), opts=getOptions();
      if(!LA.isSymmetric(A))throw new Error("La matriz A debe ser simétrica (A = Aᵀ).");
      LA.cholesky(A);
      const ref=LA.solveCholesky(A,b);
      const gc=LA.conjugateGradient(A,b,opts);
      const gcp=$("preconditioner").value==="jacobi"?LA.pcg(A,b,opts):LA.conjugateGradient(A,b,opts);
      const kA=LA.conditionNumber(A), kM=LA.preconditionedCondition(A);
      last={A,b,ref,gc,gcp,kA,kM,pre:$("preconditioner").value};
      $("validationMessage").textContent="Cálculo completado. La matriz cumple simetría y definición positiva.";
      $("validationMessage").className="validation-message success";
      renderResults(last); drawChart(gc.history,gcp.history);
    }catch(e){$("validationMessage").textContent=e.message;$("validationMessage").className="validation-message error";}
  }
  function renderResults(r){
    const gc=r.gc,gcp=r.gcp;
    $("summaryTable").innerHTML=[
      ["Solución (x₁, x₂, x₃, x₄)",`[${gc.x.map(v=>fmt(v,3)).join(", ")}]`,`[${gcp.x.map(v=>fmt(v,3)).join(", ")}]`,`[${r.ref.map(v=>fmt(v,3)).join(", ")}]`],
      ["Iteraciones",gc.iterations,gcp.iterations,"—"],
      ["Residuo relativo final",sci(gc.residual),sci(gcp.residual),"0 (referencia numérica)"],
      ["Número de condición κ(A)",fmt(r.kA,3),fmt(r.kA,3),fmt(r.kA,3)],
      ["κ precondicionada (Jacobi)",fmt(r.kM,3),fmt(r.kM,3), "—"],
      ["¿Convergió?",gc.converged?"Sí":"No",gcp.converged?"Sí":"No","—"]
    ].map(row=>`<tr>${row.map(v=>`<td>${v}</td>`).join("")}</tr>`).join("");
    $("solutionTable").innerHTML=r.ref.map((v,i)=>`<tr><td>${units[i]}</td><td>${fmt(v,6)}</td><td>${fmt(gcp.x[i],6)}</td></tr>`).join("");
    $("errorTable").innerHTML=r.ref.map((v,i)=>{
      const ea=Math.abs(v-gcp.x[i]), er=v===0?null:ea/Math.abs(v);
      return `<tr><td>x${i+1}</td><td>${sci(ea)}</td><td>${er===null?"No definido":sci(er)}</td><td>${er===null?"No definido":fmt(er*100,6)+"%"}</td></tr>`;
    }).join("");
    const label=current==="ideal"?"estabilidad y coordinación habitual":current==="stress"?"presión presupuestaria y mayores metas de servicio":"diferencias de escala y acoplamiento entre unidades";
    $("interpretation").innerHTML=`<b>Interpretación financiera:</b> el escenario representa ${label}. La solución de referencia y la aproximación iterativa indican la asignación ponderada que satisface el modelo. El residuo relativo GCP es <b>${sci(gcp.residual)}</b>; ${gcp.converged?"se alcanzó la tolerancia solicitada.":"no se alcanzó la tolerancia dentro del máximo de iteraciones."} Los valores son resultados de un modelo sintético, no recomendaciones de gasto real.`;
  }
  function drawChart(a,b){
    const c=$("convergenceChart"),ctx=c.getContext("2d"),w=c.width,h=c.height,pad={l:58,r:20,t:20,b:35};
    ctx.clearRect(0,0,w,h);ctx.fillStyle="#fff";ctx.fillRect(0,0,w,h);
    const all=a.concat(b).filter(v=>v>0), min=Math.max(1e-16,Math.min(...all)), max=Math.max(1,...all);
    const xMax=Math.max(a.length,b.length,2)-1, y=val=>pad.t+(Math.log10(max)-Math.log10(Math.max(min,val)))/(Math.log10(max)-Math.log10(min)||1)*(h-pad.t-pad.b), x=i=>pad.l+i/(xMax||1)*(w-pad.l-pad.r);
    ctx.strokeStyle="#e8edf3";ctx.lineWidth=1;ctx.fillStyle="#788597";ctx.font="12px sans-serif";
    for(let k=0;k<=4;k++){const yy=pad.t+k*(h-pad.t-pad.b)/4;ctx.beginPath();ctx.moveTo(pad.l,yy);ctx.lineTo(w-pad.r,yy);ctx.stroke();ctx.fillText((Math.pow(10,Math.log10(max)-k*(Math.log10(max)-Math.log10(min))/4)).toExponential(0),4,yy+4);}
    function line(arr,color){ctx.strokeStyle=color;ctx.lineWidth=3;ctx.beginPath();arr.forEach((v,i)=>{const xx=x(i),yy=y(v);if(i===0)ctx.moveTo(xx,yy);else ctx.lineTo(xx,yy);});ctx.stroke();}
    line(a,"#8b9aaa");line(b,"#4d65e8");
    ctx.fillStyle="#788597";ctx.fillText("Iteraciones",w/2-25,h-8);ctx.fillText("0",pad.l-4,h-12);ctx.fillText(String(xMax),w-pad.r-8,h-12);
  }
  function selectScenario(name){
    current=name;document.querySelectorAll(".scenario-btn").forEach(b=>b.classList.toggle("active",b.dataset.scenario===name));
    renderEditors();$("validationMessage").textContent="";run();
  }
  document.querySelectorAll(".scenario-btn").forEach(btn=>btn.addEventListener("click",()=>selectScenario(btn.dataset.scenario)));
  $("runButton").addEventListener("click",run);
  $("resetButton").addEventListener("click",()=>{renderEditors();run();});
  $("matrixEditor").addEventListener("input",updateDiagnostics);
  $("matrixEditor").addEventListener("change",run);
  $("vectorEditor").addEventListener("change",run);
  ["tolerance","maxIterations","initialVector","preconditioner"].forEach(id=>$(id).addEventListener("change",run));

  $("runQa").addEventListener("click",()=>{
    const tests=[];
    const add=(name,input,expected,obtained,ok)=>tests.push({name,input,expected,obtained,ok});
    try{
      for(const key of ["ideal","stress","ill"]){const s=scenarios[key];LA.cholesky(s.A);LA.solveCholesky(s.A,s.b);add("Carga de escenario",key,"Sistema simétrico SPD resuelto","Sistema válido; Cholesky completado",true);}
      const invalid=[[1,2,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]];
      let caught=false;try{if(!LA.isSymmetric(invalid))throw new Error("No simétrica");LA.cholesky(invalid);}catch(e){caught=true;}
      add("Validación de matriz no válida","Matriz no simétrica","Mostrar error sin detener la aplicación",caught?"Error detectado correctamente":"No se detectó el error",caught);
      const modified=LA.clone(scenarios.ideal.A);modified[0][0]+=0.25;
      const before=LA.solveCholesky(scenarios.ideal.A,scenarios.ideal.b),after=LA.solveCholesky(modified,scenarios.ideal.b);
      const changed=normDiff(before,after)>1e-10;
      add("Edición dinámica","a₁₁ incrementado en 0.25","La solución debe recalcularse",changed?"La solución cambió al modificar A":"La solución no cambió",changed);
      const s=scenarios.ill,opts={tol:1e-6,maxIter:100,x0:[0,0,0,0]};
      const gc=LA.conjugateGradient(s.A,s.b,opts),gcp=LA.pcg(s.A,s.b,opts);
      const better=gcp.converged&&(!gc.converged||gcp.iterations<=gc.iterations);
      add("Comparación de precondicionamiento","Escenario mal condicionado; Jacobi","GCP converge en menor o igual número de iteraciones que GC",`GC: ${gc.iterations} (${gc.converged?"convergió":"no convergió"}); GCP: ${gcp.iterations} (${gcp.converged?"convergió":"no convergió"})`,better);
    }catch(e){add("Ejecución de QA","Prueba interna","Completar ejecución",e.message,false);}
    $("qaTable").innerHTML=tests.map(t=>`<tr><td>${t.name}</td><td>${t.input}</td><td>${t.expected}</td><td>${t.obtained}</td><td><span class="${t.ok?"qa-pass":"qa-fail"}">${t.ok?"APROBADA":"REVISAR"}</span></td></tr>`).join("");
    const passed=tests.filter(t=>t.ok).length;
    $("qaStatus").textContent=`${passed} de ${tests.length} pruebas aprobadas`;
  });
  function normDiff(a,b){return Math.sqrt(a.reduce((s,v,i)=>s+(v-b[i])**2,0));}
  renderEditors();run();
})();