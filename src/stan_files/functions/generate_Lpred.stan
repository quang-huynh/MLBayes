
  vector generate_Lpred(int nbreaks, int count, vector Lobs, vector ss, real Lc, real Linf, real K, vector Z, vector yearZ) {
    matrix[nbreaks,count] dy;
    matrix[nbreaks+1,count] a;
    matrix[nbreaks+1,count] s;
    matrix[nbreaks+1,count] r;
  
    vector[count] denom;
    vector[count] numsum;
    vector[count] num;
  
    vector[count] Lpred;
  
    for (i in 1:nbreaks) {
      for (m in 1:count) {
	    dy[i, m] = yearZ[i] >= m ? 0 : m-yearZ[i];
	  }
    }
	
	if (nbreaks>1) {
      for (i in 1:(nbreaks-1)) {
        for (m in 1:count) dy[i,m] = dy[i,m] - dy[i+1,m];
      }
    }
    
    for(m in 1:count) {  
      denom[m] = 0.;
      numsum[m] = 0.;
    
      for(i in 1:(nbreaks+1)) {
        a[i, m] = 1;
        r[i, m] = 1.;
        if (i < nbreaks+1) s[i, m] = 1. - exp(-(Z[nbreaks + 2 - i]+K) * dy[nbreaks + 1 - i, m]);
        if (i == nbreaks+1) s[i, m] = 1.;
      
        if (i > 1) {
          for (j in 1:(i-1)) {
            a[i,m] = a[i,m] * exp(-Z[nbreaks + 2 - j] * dy[nbreaks + 1 - j, m]);
            r[i,m] = r[i,m] * exp(-(Z[nbreaks + 2 - j] + K) * dy[nbreaks + 1 - j, m]);
          }
        }
      
        if(i<nbreaks+1) denom[m] = denom[m] + a[i, m] * (1. - exp(-Z[nbreaks+2-i] * dy[nbreaks+1-i,m]))/Z[nbreaks+2-i];
        if(i==nbreaks+1) denom[m] = denom[m] + a[i, m]/Z[nbreaks + 2 - i];
        numsum[m] = numsum[m] + r[i, m] * s[i, m] / (Z[nbreaks + 2 - i] + K);
      }
    
      num[m] = Linf * (denom[m] - (1. - Lc/Linf)*numsum[m]);
      Lpred[m] = num[m]/denom[m];
    }
    return Lpred;   
  }

  