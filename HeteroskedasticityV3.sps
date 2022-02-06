* Encoding: UTF-8.

/*Breusch-pagan test for heteroskedasticity. 
/*macro created by Ahmad Daryanto (a.daryanto@lancaster.ac.uk)*/. 

set printback=off mprint=off errors=off .
preserve .
DEFINE BPK (iv = !charend('/')
/dv = !charend('/')
/robse = !charend('/')!default(3)).


MATRIX.
get mat/variables=!dv !iv  /names=nms /MISSING=OMIT.
compute n=nrow(mat).
*dv in original metrix.
compute Y=mat(:,1).
*===============================================.
*            OLS Regression of Raw Data .
*===============================================.
compute n=nrow(mat).
compute ones=make(n,1,1).
compute Y=mat(:,1).
compute X={ones,mat(:,2:ncol(mat))}.
compute b=(inv(sscp(X)))*t(X)*Y.
compute k=ncol(X).
*===computing standard error of b, t value and p-value of OLS ==.
compute e=Y-X*b.
compute e2=e(:,1)&*e(:,1).
compute sser=csum(e2).
compute mse=(1/(n-k))*sser.
compute vb=mse*inv(sscp(X)).
compute sb=sqrt(diag(vb)).
compute tb=b/sb.
compute dff=n-k.
compute F=tb&*tb.
compute pF=1-fcdf(F,1,dff).
compute pF=1-fcdf(F,1,dff).
     *--95% CI--.

   compute tcrit = IDF.T(0.975, dff).     
   compute LB=b - tcrit*sb.
   compute UB=b+ tcrit*sb.
   compute olsout={b,sb, tb,pF,LB,UB}.

*/predicted value.
compute yfitted = X*b.
compute resid = Y - yfitted.
compute dataplot={yfitted,resid}.
save dataplot/outfile=* /var  yfitted resid.    

*print yfitted/format f9.2.
*print resid /format f9.2.

*for output with robust std error HC0.  
do if (!robse=0).
    compute  vbh=inv(sscp(X))*t(X)*mdiag(e2)*X*inv(sscp(X)).
end if.
* HC1.
do if (!robse=1).
compute  vbh=inv(sscp(X))*t(X)*mdiag(e2)*X*inv(sscp(X)).
    compute vbh=vbh*N/(N-k).
end if.
*HC2.
do if (!robse=2).
compute hat=X*inv(sscp(X))* t(X).
compute dhat=e2&/(ones-diag(hat)).
compute vbh=inv(sscp(X))*t(X)*mdiag(dhat)*X*inv(sscp(X)).
end if.
*HC3.
do if (!robse=3).
compute hat=X*inv(sscp(X))* t(X).
compute hat2=(ones-diag(hat))&*(ones-diag(hat)).
compute dhat=e2&/hat2.
compute vbh=inv(sscp(X))*t(X)*mdiag(dhat)*X*inv(sscp(X)).
end if.
*HC4.
do if (!robse=4).
compute hat=X*inv(sscp(X))* t(X).
compute fours=make(n,1,4).
compute mh={fours,n*diag(hat)/k}.
compute dummy=rmin(mh).
compute hat2=(ones-diag(hat))&**dummy.
compute dhat=e2&/hat2.
compute vbh=inv(sscp(X))*t(X)*mdiag(dhat)*X*inv(sscp(X)).
end if.
    compute sbh=sqrt(diag(vbh)).
    compute tbh=b/sbh.
    compute dff=n-k.
    compute Fh=tbh&*tbh.
    compute pFh=1-fcdf(Fh,1,dff).
    compute pF=1-fcdf(Fh,1,dff).
     *--95% CI--.
     
   compute LBh=b - tcrit*sbh.
   compute UBh=b + tcrit*sbh.
   compute olsouth={b,sbh, tbh,pFh, LBh, UBh}.
   
*end of calculation.
print/title=" written by Ahmad Daryanto".
compute temp=t(nms(:,1)).
print/title="Original Regression model:".
print temp/title="Dependent variable"/format=A8.
*===Preparing input ANOVA table.
*computing mean square regression.
compute meanY=ones*t(csum(Y)/n).
compute e_reg=X*b-meanY.
compute ssreg=csum(sscp(e_reg)).
compute sumsq=T({ssreg,sser}).
compute dfa=T({k-1,n-k}).
compute mse_a=sumsq/dfa.
compute Fval=(ssreg/(k-1))/(sser/(n-k)).
Compute pF_a=1-fcdf(Fval,k-1,n-k).
compute F_a=T({Fval,-999}).
compute pFa=T({pF_a,-999}).
*--computing R-square.
Compute total=sser+ssreg.
Compute Rsq=ssreg/total.
print Rsq/title="R-square"/format=F9.3.
*--OLS output.
compute nmvars = t(nms(1,2:ncol(mat))).
compute nmvars = {"constant"; nmvars; "interact"}.
compute cnms={"b","se", "t", "sig", "95%LB", "95%UB"}.
print olsout/title ="OLS outputs"/rnames=nmvars/cnames=cnms/format=F9.4.


*--OLS output associated with robust standard errors.
compute nmvars = t(nms(1,2:ncol(mat))).
compute nmvars = {"constant"; nmvars; "interact"}.
print olsouth/title ="OLS outputs with heteroskedasticity-robust standard "+
    "errors:"/rnames=nmvars/cnames=cnms/format=F9.4.
do if (!robse=0).
print/title="* Note: standard error is HC0 variant (Eicker-Huber–White standard errors), not "+
    "recommended for sample sizes < 250 (Long and Ervin, 2000)".
end if.
do if (!robse=1).
print/title="* Note: standard error is HC1 variant".
end if.
do if (!robse=2).
print/title="* Note: standard error is HC2 variant".
end if.
do if (!robse=3).
print/title="* Note: standard error is HC3 variant".
end if.
do if (!robse=4).
print/title="* Note: standard error is HC4 variant".
end if.
*--ANOVA table.
print {sumsq,dfa, mse_a,F_a,pFa} /space=3
  /title '------- ANOVA TABLE  --------'
  /clabel "SS" "df" "MS" "F" "Sig"
  /rlabel "Model" "Residual"
  /format f10.4 .
/*=========================================.
/*Breusch-Pagan test for heteroskedasticity
/*=========================================.
compute var_e=sscp(e)/n.
*residuals are scaled.
compute g=e2/var_e.
compute bp=(inv(sscp(X)))*t(X)*g.
compute ep=g-X*bp.
compute e2p=ep(:,1)&*ep(:,1).
compute sserp=csum(e2p).
compute msep=(1/(n-k))*sserp.
compute vbp=msep*inv(sscp(X)).
compute sbp=sqrt(diag(vbp)).
compute tbp=bp/sbp.
compute dff=n-k.
compute Fp=tbp&*tbp.
compute pFp=1-fcdf(Fp,1,dff).
     *--95% CI--.
     
   compute LB=bp  - tcrit * sbp.
   compute UB=bp + tcrit * sbp.
   compute olsout={b,sb, tb,pF,LB,UB}.

compute olsout={bp,sbp, tbp,pFp, LB, UB}.
print/title="============================================".
print/title="Breusch-Pagan and Koenker test".
print/title="============================================".
print/title="The tests use the scaled residuals from the original OLS above with no adjustment to "+
    ""+    
    "standard errors.".
print olsout/title ="OLS outputs"/rnames=nmvars/cnames=cnms/format=F9.4.
*--Computing LM statistics .
compute meanY=ones*t(csum(g)/n).
compute e_regp=X*bp-meanY.
compute ssregp=csum(sscp(e_regp)).
Compute total=sserp+ssregp.
Compute Rsqp=ssregp/total.
print Rsqp/title="R-square"/format=F9.3.
compute F=(ssregp/(k-1))/(sserp/(n-k)).
Compute pF=1-fcdf(Fval,k-1,n-k).
*--ANOVA table.
compute F_a=T({F,-999}).
compute pF_a=T({pF,-999}).
compute sumsq=T({ssregp,sserp}).
compute msep=sumsq/dfa.
print {sumsq,dfa, msep,F_a,pF_a} /space=3
  /title '------- ANOVA TABLE  --------'
  /clabel "SS" "df" "MS" "F" "Sig"
  /rlabel "Model" "Residual"
  /format f10.4 .
/* test statisticsby Breusch-Pagan.
compute np=ncol(mat)-1.
Compute LMb=0.5*ssregp.
compute sigb=1-chicdf(LMb,np).
/* test statisticsby Koenker.
Compute LMk=n*Rsqp.
compute sigk=1-chicdf(LMk,np).
compute LM=T({LMb,LMk}).
compute sig=T({sigb,sigk}).
print{LM,sig}
  /title '------- Breusch-Pagan and Koenker test statistics and sig-values --------'
  /clabel "LM"  "Sig"
  /rlabel "BP" "Koenker"
  /format f10.4 .
print/title="Null hypothesis: heteroskedasticity not present (homoskedasticity).".
print/title="If sig-value less than 0.05, reject the null hypothesis.". 
print/title="Note: Breusch-Pagan test is a large sample test and assumes the residuals to be "+
    "normally distributed.".
    
END MATRIX.

GRAPH
  /SCATTERPLOT(BIVAR)=yfitted WITH resid
  /MISSING=LISTWISE.


!ENDDEFINE.
restore
