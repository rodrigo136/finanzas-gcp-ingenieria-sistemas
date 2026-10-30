/* Operaciones matriciales y métodos iterativos. JavaScript puro, sin dependencias. */
(function () {
  const clone = A => A.map(row => row.slice());
  const dot = (a,b) => a.reduce((s,v,i)=>s+v*b[i],0);
  const norm = a => Math.sqrt(dot(a,a));
  const matVec = (A,x) => A.map(row=>dot(row,x));
  const sub = (a,b) => a.map((v,i)=>v-b[i]);
  const addScaled = (a,b,s) => a.map((v,i)=>v+s*b[i]);
  const scale = (a,s) => a.map(v=>v*s);
  const transpose = A => A[0].map((_,j)=>A.map(row=>row[j]));

  function isSymmetric(A, tol=1e-10) {
    return A.every((row,i)=>row.every((v,j)=>Math.abs(v-A[j][i]) <= tol*Math.max(1,Math.abs(v),Math.abs(A[j][i]))));
  }

  // Cholesky comprueba definición positiva y permite resolver Ax=b.
  function cholesky(A) {
    const n=A.length, L=Array.from({length:n},()=>Array(n).fill(0));
    for(let i=0;i<n;i++) for(let j=0;j<=i;j++) {
      let sum=A[i][j];
      for(let k=0;k<j;k++) sum-=L[i][k]*L[j][k];
      if(i===j) {
        if(!(sum>1e-12) || !Number.isFinite(sum)) throw new Error("La matriz no es definida positiva (Cholesky encontró un pivote no positivo).");
        L[i][j]=Math.sqrt(sum);
      } else L[i][j]=sum/L[j][j];
    }
    return L;
  }

  function solveCholesky(A,b) {
    const L=cholesky(A), n=b.length, y=Array(n).fill(0), x=Array(n).fill(0);
    for(let i=0;i<n;i++){let s=b[i];for(let j=0;j<i;j++)s-=L[i][j]*y[j];y[i]=s/L[i][i];}
    for(let i=n-1;i>=0;i--){let s=y[i];for(let j=i+1;j<n;j++)s-=L[j][i]*x[j];x[i]=s/L[i][i];}
    return x;
  }

  // Autovalores de una matriz simétrica mediante rotaciones de Jacobi.
  function eigenvaluesSymmetric(A) {
    const B=clone(A), n=B.length, max=100*n*n;
    for(let iter=0;iter<max;iter++){
      let p=0,q=1,maxVal=0;
      for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)if(Math.abs(B[i][j])>maxVal){maxVal=Math.abs(B[i][j]);p=i;q=j;}
      if(maxVal<1e-12)break;
      const tau=(B[q][q]-B[p][p])/(2*B[p][q]);
      const t=(tau>=0?1:-1)/(Math.abs(tau)+Math.sqrt(1+tau*tau));
      const c=1/Math.sqrt(1+t*t), s=t*c;
      const app=B[p][p], aqq=B[q][q], apq=B[p][q];
      B[p][p]=app-t*apq; B[q][q]=aqq+t*apq; B[p][q]=B[q][p]=0;
      for(let k=0;k<n;k++)if(k!==p&&k!==q){const bkp=B[k][p],bkq=B[k][q];B[k][p]=B[p][k]=c*bkp-s*bkq;B[k][q]=B[q][k]=s*bkp+c*bkq;}
    }
    return B.map((row,i)=>row[i]).sort((a,b)=>a-b);
  }

  function conditionNumber(A) {
    const ev=eigenvaluesSymmetric(A), min=ev[0], max=ev[ev.length-1];
    return min>0 ? max/min : Infinity;
  }

  function residual(A,x,b){return sub(b,matVec(A,x));}
  function relativeResidual(A,x,b){const nb=norm(b);return nb===0?norm(residual(A,x,b)):norm(residual(A,x,b))/nb;}

  function conjugateGradient(A,b,opts={}) {
    const n=b.length, tol=opts.tol??1e-6, maxIter=opts.maxIter??100;
    let x=(opts.x0??Array(n).fill(0)).slice();
    let r=residual(A,x,b), p=r.slice(), rr=dot(r,r), nb=norm(b)||1;
    const history=[norm(r)/nb]; let converged=history[0]<tol, iterations=0;
    while(!converged && iterations<maxIter){
      const Ap=matVec(A,p), denom=dot(p,Ap);
      if(!(denom>0) || !Number.isFinite(denom)) break;
      const alpha=rr/denom;
      x=addScaled(x,p,alpha); r=addScaled(r,Ap,-alpha);
      const rrNew=dot(r,r); iterations++;
      history.push(Math.sqrt(rrNew)/nb);
      converged=history[history.length-1]<tol;
      if(converged){rr=rrNew;break;}
      const beta=rrNew/rr; p=addScaled(r,p,beta); rr=rrNew;
    }
    return {x,iterations,history,residual:relativeResidual(A,x,b),converged};
  }

  function pcg(A,b,opts={}) {
    const n=b.length, tol=opts.tol??1e-6, maxIter=opts.maxIter??100;
    const x0=(opts.x0??Array(n).fill(0)).slice();
    const invDiag=A.map((row,i)=>1/row[i]);
    if(invDiag.some(v=>!Number.isFinite(v)||v<=0)) throw new Error("Jacobi requiere elementos diagonales positivos.");
    let x=x0, r=residual(A,x,b), z=r.map((v,i)=>v*invDiag[i]), p=z.slice(), rz=dot(r,z), nb=norm(b)||1;
    const history=[norm(r)/nb]; let converged=history[0]<tol, iterations=0;
    while(!converged && iterations<maxIter){
      const Ap=matVec(A,p), denom=dot(p,Ap);
      if(!(denom>0) || !Number.isFinite(denom)) break;
      const alpha=rz/denom;
      x=addScaled(x,p,alpha); r=addScaled(r,Ap,-alpha);
      const rel=norm(r)/nb; iterations++; history.push(rel); converged=rel<tol;
      if(converged)break;
      z=r.map((v,i)=>v*invDiag[i]);
      const rzNew=dot(r,z), beta=rzNew/rz;
      p=addScaled(z,p,beta); rz=rzNew;
    }
    return {x,iterations,history,residual:relativeResidual(A,x,b),converged};
  }

  function preconditionedCondition(A) {
    // M=diag(A); M^(-1/2) A M^(-1/2) es simétrica y tiene
    // los mismos autovalores que M^(-1)A.
    const d=A.map((row,i)=>row[i]);
    const B=A.map((row,i)=>row.map((v,j)=>v/Math.sqrt(d[i]*d[j])));
    return conditionNumber(B);
  }

  window.LinearAlgebra={clone,dot,norm,matVec,sub,isSymmetric,cholesky,solveCholesky,eigenvaluesSymmetric,conditionNumber,relativeResidual,conjugateGradient,pcg,preconditionedCondition};
})();