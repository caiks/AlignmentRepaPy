from AlignmentRepa import *
from AlignmentRandomRepa import *

# parametersSystemsBuilderTupleNoSumlayerMultiEffectiveRepa_ui :: 
#   Integer -> Integer -> Integer -> Integer -> System -> Set.Set Variable -> Fud -> 
#   HistoryRepa -> HistogramRepaRed -> HistoryRepa -> HistogramRepaRed ->   
#   ([((Set.Set Variable, (HistogramRepaVec, HistogramRepaVec, UV.Vector Double)),Double)],Integer)

def parametersSystemsBuilderTupleNoSumlayerMultiEffectiveRepa_ui(xmax,omax,bmax,mmax,uu,vv,ff,hh,hhx,hhrr,hhrrx):
    def sgl(x):
        return sset([x])
    def top(amax,ll):
        ll1 = ll.copy()
        ll1.sort(key = lambda x: x[0], reverse = True)
        return ll1[:amax]
    def topd(amax,ll):
        return [b for (a,b) in top(amax,ll)]
    size = historyRepasSize
    vrrrrv = vectorHistogramRepasHistogramRepaVec_u
    xind = histogramRepaRedsIndependent
    def xred(hhx,vv):
        return setVarsHistogramRepaRedsRed(vv,hhx)
    reduce = setVarsHistoryRepasReduce
    cross = parametersSetVarsHistoryRepasSetSetVarsAlignedTop_u
    append = parametersSetVarsSetSetVarsHistoryRepasSetSetVarsAlignedTop_u
    rrvffv = histogramRepaVecsFaclnsRepaVecs
    rrvsum = histogramRepaVecsSum
    fder = fudsDerived
    fvars = fudsVars
    def meff(hhx,vv):
        (_,mvv,_,lrr) = hhx
        return sset([v for v in vv if len([a for a in np.nditer(lrr[mvv[v]]) if a != 0]) > 1])
    z = size(hh)
    zrr = size(hhrr)
    f = float(z/zrr)
    vv1 = meff(hhx,vv)
    if len(vv1) < 2:
        return ([],0)
    def buildb(ww,qq,nn,sn):
        if len(qq) == 0:
            return (nn,sn)
        qq1 = [b for (a,b) in qq]
        (mm,sm) = append(xmax,omax,ww,qq1,hh,hhx,hhrr,hhrrx)
        if len(mm) > 0:
            return buildb(ww,mm,nn+mm,sn+sm) 
        return (nn,sn)
    def res(xx):
        x2 = []
        for jj in xx:
            bb = reduce(1,jj,hh)
            bbrr = reduce(f,jj,hhrr)
            bbx = xind(z,xred(hhx,jj))
            bbrrx = xind(z,xred(hhrrx,jj))
            bbv = vrrrrv(z,[bb,bbx,bbrr,bbrrx])
            ffv = rrvffv(bbv)
            ssv = rrvsum(ffv)
            [a1,_,b1,_] = ssv
            x2.append(((jj,(bbv,ffv,ssv)),a1-b1))
        return x2
    if len(ff) == 0:
        (xc,sc) = cross(xmax,omax,vv1,hh,hhx,hhrr,hhrrx)
        (x0,s0) = buildb(vv1,xc,xc,sc)
        return (res(topd(bmax//mmax,x0)),s0)
    yy = fvars(ff)|vv1
    vdd = [sgl(w) for w in fder(ff)]
    (xa,sa) = append(xmax,omax,yy,vdd,hh,hhx,hhrr,hhrrx)
    (x1,s1) = buildb(yy,xa,xa,sa)
    return (res(topd(bmax//mmax,x1)),s1)

# parametersSystemsBuilderTupleNoSumlayerMultiEffectiveRepa_ui_1 :: 
#   Integer -> Integer -> Integer -> Integer -> System -> Set.Set Variable -> Fud -> 
#   HistoryRepa -> HistogramRepaRed -> HistoryRepa -> HistogramRepaRed ->   
#   ([((Set.Set Variable, (HistogramRepaVec, HistogramRepaVec, UV.Vector Double)),Double)],Integer)

def parametersSystemsBuilderTupleNoSumlayerMultiEffectiveRepa_ui_1(xmax,omax,bmax,mmax,uu,vv,ff,hh,hhx,hhrr,hhrrx):
    def sgl(x):
        return sset([x])
    def top(amax,ll):
        ll1 = ll.copy()
        ll1.sort(key = lambda x: x[0], reverse = True)
        return ll1[:amax]
    def topd(amax,ll):
        return [b for (a,b) in top(amax,ll)]
    vol = systemsSetVarsVolume_u
    size = historyRepasSize
    vrrrrv = vectorHistogramRepasHistogramRepaVec_u
    hvempty = histogramRepaVecEmpty
    xind = histogramRepaRedsIndependent
    def xred(hhx,vv):
        return setVarsHistogramRepaRedsRed(vv,hhx)
    reduce = setVarsHistoryRepasReduce
    rrvffv = histogramRepaVecsFaclnsRepaVecs
    rrvsum = histogramRepaVecsSum
    fder = fudsDerived
    fvars = fudsVars
    def meff(hhx,vv):
        (_,mvv,_,lrr) = hhx
        return sset([v for v in vv if len([a for a in np.nditer(lrr[mvv[v]]) if a != 0]) > 1])
    def init(vv):
        return [((0,0,0),((sgl(w),(hvempty(),hvempty(),[])),0)) for w in vv]
    def final(nn):
        return [(a,((kk,xx),y)) for (a,((kk,xx),y)) in nn if len(kk) > 1]
    z = size(hh)
    zrr = size(hhrr)
    f = float(z/zrr)
    def buildb(ww,qq,nn,s2):
        pp = sset([kk|sgl(w) for (_,((kk,_),_)) in qq for w in (ww-kk)])
        x2 = []
        for jj in pp:
            u = vol(uu,jj)
            if u <= xmax:
                bb = reduce(1,jj,hh)
                bbrr = reduce(f,jj,hhrr)
                bbx = xind(z,xred(hhx,jj))
                bbrrx = xind(z,xred(hhrrx,jj))
                bbv = vrrrrv(z,[bb,bbx,bbrr,bbrrx])
                ffv = rrvffv(bbv)
                ssv = rrvsum(ffv)
                [a1,a2,b1,b2] = ssv
                x2.append(((a1-a2-b1+b2,-b1+b2,-u),((jj,(bbv,ffv,ssv)),a1-b1)))
        mm = top(omax,x2)
        if len(mm) > 0:
            return buildb(ww,mm,nn+mm,s2+len(x2)) 
        return (final(nn),s2)
    vv1 = meff(hhx,vv)
    if len(vv1) == 0:
        return ([],0)
    if len(ff) == 0:
        (x0,s0) = buildb(vv1,init(vv1),[],0)
        return (topd(bmax//mmax,x0),s0)
    (x1,s1) = buildb(fvars(ff)|vv1,init(fder(ff)),[],0)
    return (topd(bmax//mmax,x1),s1)

# parametersSystemsBuilderTupleNoSumlayerMultiEffectiveRepa_ui_2 :: 
#   Integer -> Integer -> Integer -> Integer -> System -> Set.Set Variable -> Fud -> 
#   HistoryRepa -> HistogramRepaRed -> HistoryRepa -> HistogramRepaRed ->   
#   ([((Set.Set Variable, (HistogramRepaVec, HistogramRepaVec, UV.Vector Double)),Double)],Integer)

def parametersSystemsBuilderTupleNoSumlayerMultiEffectiveRepa_ui_2(xmax,omax,bmax,mmax,uu,vv,ff,hh,hhx,hhrr,hhrrx):
    def sgl(x):
        return sset([x])
    def top(amax,ll):
        ll1 = ll.copy()
        ll1.sort(key = lambda x: x[0], reverse = True)
        return ll1[:amax]
    def topd(amax,ll):
        return [b for (a,b) in top(amax,ll)]
    vol = systemsSetVarsVolume_u
    size = historyRepasSize
    vrrrrv = vectorHistogramRepasHistogramRepaVec_u
    hvempty = histogramRepaVecEmpty
    xind = histogramRepaRedsIndependent
    def xred(hhx,vv):
        return setVarsHistogramRepaRedsRed(vv,hhx)
    reduce = setVarsHistoryRepasReduce
    rrvffv = histogramRepaVecsFaclnsRepaVecs
    rrvsum = histogramRepaVecsSum
    fder = fudsDerived
    fvars = fudsVars
    def meff(hhx,vv):
        (_,mvv,_,lrr) = hhx
        return sset([v for v in vv if len([a for a in np.nditer(lrr[mvv[v]]) if a != 0]) > 1])
    def init(vv):
        return [((0,0,0),((sgl(w),(hvempty(),hvempty(),[])),0)) for w in vv]
    def final(nn):
        return [(a,((kk,xx),y)) for (a,((kk,xx),y)) in nn if len(kk) > 1]
    z = size(hh)
    zrr = size(hhrr)
    f = float(z/zrr)
    def buildb(ww,qq,nn,s2):
        pp = sset([kk|sgl(w) for (_,((kk,_),_)) in qq for w in (ww-kk)])
        x2 = []
        for jj in pp:
            u = vol(uu,jj)
            if u <= xmax:
                bb = reduce(1,jj,hh)
                bbrr = reduce(f,jj,hhrr)
                bbx = xind(z,xred(hhx,jj))
                bbrrx = xind(z,xred(hhrrx,jj))
                bbv = vrrrrv(z,[bb,bbx,bbrr,bbrrx])
                ffv = rrvffv(bbv)
                ssv = rrvsum(ffv)
                [a1,a2,b1,b2] = ssv
                x2.append(((a1-a2-b1+b2,-b1+b2,-u),((jj,(bbv,ffv,ssv)),a1-b1)))
        mm = top(omax,x2)
        if len(mm) > 0:
            s3 = s2+len(x2)
            x2 = None
            return buildb(ww,mm,nn+mm,s3) 
        return (final(nn),s2)
    vv1 = meff(hhx,vv)
    if len(vv1) == 0:
        return ([],0)
    if len(ff) == 0:
        (x0,s0) = buildb(vv1,init(vv1),[],0)
        return (topd(bmax//mmax,x0),s0)
    (x1,s1) = buildb(fvars(ff)|vv1,init(fder(ff)),[],0)
    return (topd(bmax//mmax,x1),s1)

if not AlignmentForeignPy_ok:
    parametersSystemsBuilderTupleNoSumlayerMultiEffectiveRepa_ui = parametersSystemsBuilderTupleNoSumlayerMultiEffectiveRepa_ui_2

# parametersSystemsBuilderTupleNoSumlayerMultiEffectiveRepa_u :: 
#   Integer -> Integer -> Integer -> Integer -> System -> Set.Set Variable -> Fud -> 
#   HistoryRepa -> HistogramRepaRed -> HistoryRepa -> HistogramRepaRed ->   
#   [((Set.Set Variable, (HistogramRepaVec, HistogramRepaVec, UV.Vector Double)),Double)]

def parametersSystemsBuilderTupleNoSumlayerMultiEffectiveRepa_u(xmax,omax,bmax,mmax,uu,vv,ff,hh,hhx,hhrr,hhrrx):
    return parametersSystemsBuilderTupleNoSumlayerMultiEffectiveRepa_ui(xmax,omax,bmax,mmax,uu,vv,ff,hh,hhx,hhrr,hhrrx)[0]

# parametersSystemsPartitionerMaxRollByMRepa_ui :: 
#   Integer -> Integer -> Integer -> System -> Set.Set Variable -> 
#   (HistogramRepaVec, HistogramRepaVec, UV.Vector Double) -> Double ->
#   ([(Set.Set (Set.Set Variable),(HistogramRepaVec,HistogramRepaVec))],Integer)

def parametersSystemsPartitionerMaxRollByMRepa_ui(mmax,umax,pmax,uu,kk,hr,y1):
    rrvvrr = histogramRepaVecsArray
    pprr = setSetVarsHistogramRepaVecsPartitionVec_u
    ppxx = setSetVarsHistogramRepaPairPartitionIndependentPair_u
    rrvqqy = parametersHistogramRepaVecsSetTuplePartitionTopByM_u
    (rrv,_,_) = hr
    (vbb,mbb,z,sbb,rbb) = rrv
    [bb,_,bbrr,_] = rbb
    nnv = (vbb,mbb,z,sbb,[bb,bbrr])
    n = len(sbb)
    (mm2,q) = rrvqqy(mmax,umax,pmax,rrv,y1)
    mm3 = []
    for yy in mm2:
        ccv = pprr(yy,nnv)
        (vcc,mcc,_,scc,rcc) = ccv
        [cc,ccrr] = rcc
        [ccx,ccrrx] = rrvvrr(ppxx(yy,nnv))
        ff = sp.gammaln(cc + 1)
        ffrr = sp.gammaln(ccrr + 1)
        ffx = sp.gammaln(ccx + 1)
        ffrrx = sp.gammaln(ccrrx + 1)
        mm3.append((yy,((vcc,mcc,z,scc,[cc,ccx,ccrr,ccrrx]),(vcc,mcc,1,scc,[ff,ffx,ffrr,ffrrx]))))
    return (mm3,q)

# parametersSystemsPartitionerMaxRollByMRepa_ui_1 :: 
#   Integer -> Integer -> Integer -> System -> Set.Set Variable -> 
#   (HistogramRepaVec, HistogramRepaVec, UV.Vector Double) -> Double ->
#   ([(Set.Set (Set.Set Variable),(HistogramRepaVec,HistogramRepaVec))],Integer)

def parametersSystemsPartitionerMaxRollByMRepa_ui_1(mmax,umax,pmax,uu,kk,bb,y1):
    def top(amax,ll):
        ll1 = ll.copy()
        ll1.sort(key = lambda x: x[0], reverse = True)
        return ll1[:amax]
    def topd(amax,ll):
        return [b for (a,b) in top(amax,ll)]
    prod = listSetsProduct
    stirsll = setsSetPartitionFixed
    vol = systemsSetVarsVolume_u
    rrvvrr = histogramRepaVecsArray
    rrvffv = histogramRepaVecsFaclnsRepaVecs
    rrvsum = histogramRepaVecsSum
    pprr = setSetVarsHistogramRepaVecsPartitionVec_u
    ppxx = setSetVarsHistogramRepaVecsPartitionIndependentVec_u
    s = 0
    (rrv,_,_) = bb
    (vbb,mbb,z,sbb,rbb) = rrv
    [bb,bbx,bbrr,bbrrx] = rbb
    nnv = (vbb,mbb,z,sbb,[bb,bbrr])
    n = len(sbb)
    mmax1 = mmax
    if mmax1 >= n:
        mmax1 = n
    v = vol(uu,kk)
    pp = []
    for m in range(2,mmax1+1):
        c = v ** (1.0/m)
        mm = []
        for yy in stirsll(kk,m):
            if all([vol(uu,jj) <= umax for jj in yy]):
                ccv = ppxx(yy,nnv)
                ffv = rrvffv(ccv)
                [a2,b2] = rrvsum(ffv)
                mm.append((((y1-a2+b2)/c,b2,-m),(yy,(ccv,ffv))))
                s += 1
        pp.extend(topd(pmax,mm))
    m3 = []
    for (yy,((vcc,mcc,z,scc,rcc),(vff,mff,_,sff,rff))) in pp:
        [ccx,ccrrx] = rcc
        [cc,ccrr] = rrvvrr(pprr(yy,nnv))
        [ffx,ffrrx] = rff
        ff = sp.gammaln(cc + 1)
        ffrr = sp.gammaln(ccrr + 1)
        m3.append((yy,((vcc,mcc,z,scc,[cc,ccx,ccrr,ccrrx]),(vff,mff,1,sff,[ff,ffx,ffrr,ffrrx]))))
    return (m3,s)

if not AlignmentForeignPy_ok:
    parametersSystemsPartitionerMaxRollByMRepa_ui = parametersSystemsPartitionerMaxRollByMRepa_ui_1

# parametersSystemsPartitionerMaxRollByMRepa_u :: 
#   Integer -> Integer -> Integer -> System -> Set.Set Variable -> 
#   (HistogramRepaVec, HistogramRepaVec, UV.Vector Double) -> Double ->
#   [(Set.Set (Set.Set Variable),(HistogramRepaVec,HistogramRepaVec))]

def parametersSystemsPartitionerMaxRollByMRepa_u(mmax,umax,pmax,uu,kk,bb,y1):
    return parametersSystemsPartitionerMaxRollByMRepa_ui(mmax,umax,pmax,uu,kk,bb,y1)[0]

# parametersRollerMaximumRollExcludedSelfRepa_i :: 
#   (Set.Set (Set.Set Variable),(HistogramRepaVec,HistogramRepaVec)) -> 
#   ([(Set.Set (Set.Set Variable),V.Vector (UV.Vector Int))],Integer)

def parametersRollerMaximumRollExcludedSelfRepa_i(qq):
    (yy,(rrv,_)) = qq
    (tt,q) = histogramRepaVecsRollMax(rrv)
    for vv in tt:
        if max(vv) < len(vv) - 1:
            return ([(yy,tt)],q)
    return  ([],q)

# parametersRollerMaximumRollExcludedSelfRepa_i_1 :: 
#   (Set.Set (Set.Set Variable),(HistogramRepaVec,HistogramRepaVec)) -> 
#   ([(Set.Set (Set.Set Variable),V.Vector (UV.Vector Int))],Integer)

def parametersRollerMaximumRollExcludedSelfRepa_i_1(qq):
    def prod(ll):
        p = 1
        for x in ll:
            p *= x
        return p
    def top(ll):
        if len(ll) > 0:
            (a,b) = ll[0]
            for (a1,b1) in ll[1:]:
                if b1 > b:
                    a = a1
                    b = b1
            return [(a,b)]
        return []
    def topd(ll):
        return [a for (a,b) in top(ll)]
    def rollr(s,t,ll):
        return [r-1 if r > s else (t if r == s else r) for r in ll]
    sh = histogramRepaVecsShape
    vec = histogramRepasArray
    sing = varsHistogramRepaRedsSingle_u
    sing4 = varsHistogramRepa4VecsReduceSingle_u
    reds = histogramRepa4VecsRed_u
    faclns = histogramRepaVecsFaclnsRepaVecs
    rollv = varsHistogramRepaVecsRollVec_u
    sumv = sumsHistogramRepasRollMapPairsHistogramRepasSum_u
    copyv = varsSourcesTargetsRollsHistogramRepaVecsHistogramRepaVecRollsCopyVec_u
    (yy,(ccv,ffv)) = qq
    xxv = reds(ffv)
    a = np.sum(vec(sing(0,xxv)))
    scc = sh(ccv)
    w = prod(scc)
    m = len(scc)
    def rollb(pp,ccv,ffv,xxv,a,w,nn,y):
        scc = sh(ccv)
        mm = []
        for v in range(m):
            d = scc[v]
            if d > 2:
                w1 = w * (d-1) / d
                c1 = w1 ** (1.0/m)
                (rrv,(rs,rt)) = rollv(v,ccv)
                ggv = faclns(rrv)
                av = sumv(a,sing(v,xxv),(rs,rt),sing4(v,ggv))
                for (q,(s,t,a1)) in enumerate(zip(rs,rt,av)):
                    mm.append(((a1,w1,(v,s,t,q,rrv,ggv)),a1/c1))
        if len(mm) > 0:
            [((a1,w1,(v,s,t,q,rrv,ggv)),b1)] = top(mm)
            ccv1 = copyv(v,s,t,q,ccv,rrv)
            ffv1 = copyv(v,s,t,q,ffv,ggv)
            xxv1 = reds(ffv1)
            pp1 = pp.copy()
            pp1[v] = rollr(s,t,pp[v])
            y1 = y + len(mm)
            nn1 = nn.copy()
            nn1.append(((yy,pp1),b1))
            mm = None
            return rollb(pp1,ccv1,ffv1,xxv1,a1,w1,nn1,y1)
        return (nn,y)
    pp = [list(range(d)) for d in scc]
    (z1,y1) = rollb(pp,ccv,ffv,xxv,a,w,[],0)
    return (topd(z1),y1)

if not AlignmentForeignPy_ok:
    parametersRollerMaximumRollExcludedSelfRepa_i = parametersRollerMaximumRollExcludedSelfRepa_i_1

# parametersRollerMaximumRollExcludedSelfRepa :: 
#   (Set.Set (Set.Set Variable),(HistogramRepaVec,HistogramRepaVec)) -> 
#   [(Set.Set (Set.Set Variable),V.Vector (UV.Vector Int))]

def parametersRollerMaximumRollExcludedSelfRepa(qq):
    return parametersRollerMaximumRollExcludedSelfRepa_i(qq)[0]

# parametersSystemsBuilderDerivedVarsHighestNoSumlayerRepa_ui :: 
#   Integer -> Integer -> System -> Set.Set Variable -> Fud -> 
#   HistoryRepa -> HistogramRepaRed -> HistoryRepa -> HistogramRepaRed ->   
#   ([((Set.Set Variable, HistogramRepa, HistogramRepa), Double)],Integer)

def parametersSystemsBuilderDerivedVarsHighestNoSumlayerRepa_ui(wmax,omax,uu,vv,ff,hh,hhx,hhrr,hhrrx):
    def sgl(x):
        return sset([x])
    def top(amax,ll):
        ll1 = ll.copy()
        ll1.sort(key = lambda x: x[0], reverse = True)
        return ll1[:amax]
    def topd(amax,ll):
        return [b for (a,b) in top(amax,ll)]
    def maxd(mm):
        return [x for ((a,_,_),x) in top(1,mm)]
    vol = systemsSetVarsVolume_u
    size = historyRepasSize
    vrrrrv = vectorHistogramRepasHistogramRepaVec_u
    hvempty = histogramRepaVecEmpty
    xind = histogramRepaRedsIndependent
    def xred(hhx,vv):
        return setVarsHistogramRepaRedsRed(vv,hhx)
    reduce = setVarsHistoryRepasReduce
    append = parametersSetVarsSetSetVarsHistoryRepasSetSetVarsAlignedExcludeHiddenDenseTop_u
    rrvffv = histogramRepaVecsFaclnsRepaVecs
    rrvsum = histogramRepaVecsSum
    fder = fudsDerived
    fvars = fudsVars
    def depends(ff,w):
        return fudsSetVarsDepends(ff,sset([w]))
    yy = fvars(ff) - vv
    if len(yy) == 0:
        return ([],0)
    z = size(hh)
    zrr = size(hhrr)
    f = float(z/zrr)    
    cc = sset([(w,u) for w in (fvars(ff)-vv) for u in (fvars(depends(ff,w))-vv) if u != w])
    def buildd(ww,qq,nn,sn):
        if len(qq) == 0:
            return (nn,sn)
        qq1 = [b for (a,b) in qq]
        (mm,sm) = append(wmax,omax,cc,ww,qq1,hh,hhx,hhrr,hhrrx)
        if len(mm) > 0:
            return buildd(ww,mm,nn+mm,sn+sm) 
        return (nn,sn)
    def res(xx):
        x2 = []
        for jj in xx:
            u = vol(uu,jj)
            bb = reduce(1,jj,hh)
            bbrr = reduce(f,jj,hhrr)
            bbx = xind(z,xred(hhx,jj))
            bbrrx = xind(z,xred(hhrrx,jj))
            bbv = vrrrrv(z,[bb,bbx,bbrr,bbrrx])
            ffv = rrvffv(bbv)
            ssv = rrvsum(ffv)
            [a1,a2,b1,b2] = ssv
            m = len(jj)
            c = u ** (1.0/m)
            x2.append(((jj,bb,bbrr),(a1-a2-b1+b2)/c))
        return x2
    vdd = [sgl(w) for w in fder(ff)]
    (xa,sa) = append(wmax,omax,cc,yy,vdd,hh,hhx,hhrr,hhrrx)
    (x1,s1) = buildd(yy,xa,xa,sa)
    return (res(maxd(x1)),s1)

# parametersSystemsBuilderDerivedVarsHighestNoSumlayerRepa_ui_1 :: 
#   Integer -> Integer -> System -> Set.Set Variable -> Fud -> 
#   HistoryRepa -> HistogramRepaRed -> HistoryRepa -> HistogramRepaRed ->   
#   ([((Set.Set Variable, HistogramRepa, HistogramRepa), Double)],Integer)

def parametersSystemsBuilderDerivedVarsHighestNoSumlayerRepa_ui_1(wmax,omax,uu,vv,ff,hh,hhx,hhrr,hhrrx):
    def sgl(x):
        return sset([x])
    def top(amax,ll):
        ll1 = ll.copy()
        ll1.sort(key = lambda x: x[0], reverse = True)
        return ll1[:amax]
    def topd(amax,ll):
        return [b for (a,b) in top(amax,ll)]
    def maxfst(mm):
        return [(x,a) for ((a,_,_),x) in top(1,mm)]
    vol = systemsSetVarsVolume_u
    size = historyRepasSize
    vrrrrv = vectorHistogramRepasHistogramRepaVec_u
    hvempty = histogramRepaVecEmpty
    xind = histogramRepaRedsIndependent
    def xred(hhx,vv):
        return setVarsHistogramRepaRedsRed(vv,hhx)
    reduce = setVarsHistoryRepasReduce
    rrvffv = histogramRepaVecsFaclnsRepaVecs
    rrvsum = histogramRepaVecsSum
    fder = fudsDerived
    fvars = fudsVars
    depends = fudsSetVarsDepends
    def init(vv):
        return [((0,0,0),(sgl(w),hvempty(),hvempty())) for w in vv]
    def final(nn):
        return [(a,(kk,xx,yy)) for (a,(kk,xx,yy)) in nn if len(kk) > 1]
    z = size(hh)
    zrr = size(hhrr)
    f = float(z/zrr)
    def buildd(ww,qq,nn,s2):
        pp = sset([kk|sgl(w) for (_,(kk,_,_)) in qq for w in (ww-kk)])
        x2 = []
        for jj in pp:
            u = vol(uu,jj)
            if u <= wmax and fder(depends(ff,jj)) == jj:
                bb = reduce(1,jj,hh)
                bbrr = reduce(f,jj,hhrr)
                bbx = xind(z,xred(hhx,jj))
                bbrrx = xind(z,xred(hhrrx,jj))
                bbv = vrrrrv(z,[bb,bbx,bbrr,bbrrx])
                ffv = rrvffv(bbv)
                ssv = rrvsum(ffv)
                [a1,a2,b1,b2] = ssv
                m = len(jj)
                c = u ** (1.0/m)
                x2.append((((a1-a2-b1+b2)/c,(-b1+b2)/c,-u),(jj,bb,bbrr)))
        mm = top(omax,x2)
        if len(mm) > 0:
            s3 = s2 + len(x2)
            x2 = None
            return buildd(ww,mm,nn+mm,s3) 
        return (final(nn),s2)
    (x1,s1) = buildd(fvars(ff)-vv,init(fder(ff)),[],0)
    return (maxfst(x1),s1)

if not AlignmentForeignPy_ok:
    parametersSystemsBuilderDerivedVarsHighestNoSumlayerRepa_ui = parametersSystemsBuilderDerivedVarsHighestNoSumlayerRepa_ui_1

# parametersSystemsBuilderDerivedVarsHighestNoSumlayerRepa_u :: 
#   Integer -> Integer -> System -> Set.Set Variable -> Fud -> 
#   HistoryRepa -> HistogramRepaRed -> HistoryRepa -> HistogramRepaRed ->   
#   [((Set.Set Variable, HistogramRepa, HistogramRepa), Double)]

def parametersSystemsBuilderDerivedVarsHighestNoSumlayerRepa_u(wmax,omax,uu,vv,ff,hh,hhx,hhrr,hhrrx):
    return parametersSystemsBuilderDerivedVarsHighestNoSumlayerRepa_ui(wmax,omax,uu,vv,ff,hh,hhx,hhrr,hhrrx)[0]

# parametersSystemsLayererMaxRollByMExcludedSelfHighestRepa :: 
#   Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> 
#   System -> Set.Set Variable -> HistoryRepa -> HistogramRepaRed -> HistoryRepa -> HistogramRepaRed -> Integer ->
#   Maybe (System, Fud, [(Set.Set Variable, Double)])

def parametersSystemsLayererMaxRollByMExcludedSelfHighestRepa(wmax,lmax,xmax,omax,bmax,mmax,umax,pmax,uu,vv,xx,xxp,xxrr,xxrrp,f):
    repaRounding = 1e-6
    def sgl(x):
        return sset([x])
    def maxr(mm):
        if len(mm) > 0:
            return list(sset([b for (_,b) in mm]))[-1:][0]
        return 0
    uvars = systemsSetVar
    cart = systemsSetVarsSetStateCartesian_u
    lluu = listsSystem_u
    uunion = pairSystemsUnion
    sunion = pairStatesUnionLeft
    ssgl = stateSingleton
    llaa = listsHistogram_u
    hhvvr = historyRepasVectorVar
    apvvr = histogramRepaRedsVectorVar
    hrhx = historyRepasRed
    def unit(qq):
        return llaa([(ss,1) for ss in qq])
    tttr = systemsTransformsTransformRepa_u
    apply = historyRepasListTransformRepasApply_u
    trans = histogramsSetVarsTransform_u
    ttpp = transformsPartition
    und = transformsUnderlying
    qqff = setTransformsFud_u
    ffqq = fudsSetTransform
    def funion(ff,gg):
        return qqff(ffqq(ff) | ffqq(gg))
    def buildfftup(uu,vv,ff,hh,hhp,hhrr,hhrrp):
        return parametersSystemsBuilderTupleNoSumlayerMultiEffectiveRepa_u(xmax,omax,bmax,mmax,uu,vv,ff,hh,hhp,hhrr,hhrrp)
    def parter(uu,kk,bb,y1):
        return parametersSystemsPartitionerMaxRollByMRepa_u(mmax,umax,pmax,uu,kk,bb,y1)
    def roller(qq):
        return parametersRollerMaximumRollExcludedSelfRepa(qq)
    def buildffdervar(uu,vv,ff,xx,xxp,xxrr,xxrrp):
        return [(kk,a) for ((kk,_,_),a) in parametersSystemsBuilderDerivedVarsHighestNoSumlayerRepa_u(wmax,omax,uu,vv,ff,xx,xxp,xxrr,xxrrp)]
    def layer(vv,uu,ff,mm,xx,xxp,xxrr,xxrrp,f,l):
        if l > lmax:
            return (uu,ff,mm)
        ll0 = []
        for ((kk,bb),y1) in buildfftup(uu,vv,ff,xx,xxp,xxrr,xxrrp):
            for qq in parter(uu,kk,bb,y1):
                for (yy,pp) in roller(qq):
                    for (jj,p) in zip(yy,pp):
                        if max(p) + 1 < len(p):
                            ii = list(zip(cart(uu,jj),p))
                            ll0.append(ii)
        ll = []
        for (b,ii) in enumerate(ll0):
            w = VarPair((VarPair((VarInt(f),VarInt(l))),VarInt(b+1)))
            ww = sset([ValInt(u) for (_,u) in ii])
            tt = trans(unit([sunion(ss,ssgl(w,ValInt(u))) for (ss,u) in ii]),sgl(w))
            ll.append((tt,(w,ww)))
        ll1 = []
        for (tt,(w,ww)) in ll:
            if all([len(ww) != len(ww1) or und(tt) != und(tt1) or ttpp(tt) != ttpp(tt1) for (tt1,(w1,ww1)) in ll if w > w1]):
                ll1.append((tt,(w,ww)))
        if len(ll1) > 0:
            hh = qqff(sset([tt for (tt,_) in ll1]))
            uu1 = uunion(uu,lluu([(w,ww) for (_,(w,ww)) in ll1]))
            ffr = [tttr(uu1,tt) for (tt,_) in ll1]
            xx1 = apply(xx,ffr)
            xxp1 = hrhx(xx1)
            xxrr1 = apply(xxrr,ffr)
            xxrrp1 = hrhx(xxrr1)
            gg = funion(ff,hh)
            mm1 = buildffdervar(uu1,vv,gg,xx1,xxp1,xxrr1,xxrrp1)
            if len(mm) == 0 or maxr(mm1) > maxr(mm) + repaRounding:
                (ffr,ll0,ll,ll1) = (None,None,None,None)
                return layer(vv,uu1,gg,mm1,xx1,xxp1,xxrr1,xxrrp1,f,l+1)
        return (uu,ff,mm) 
    if wmax < 0 or lmax < 0 or xmax < 0 or omax < 0 or bmax < 0 or mmax < 1 or umax < 0 or pmax < 0:
        return None
    if not (sset(hhvvr(xx)).issubset(sset(uvars(uu))) and hhvvr(xx) == hhvvr(xxrr) and hhvvr(xx) == apvvr(xxp) and hhvvr(xx) == apvvr(xxrrp) and vv.issubset(sset(hhvvr(xx)))):
        return None
    return layer(vv,uu,fudEmpty(),[],xx,xxp,xxrr,xxrrp,f,1)

# parametersSystemsHistoryRepasDecomperMaxRollByMExcludedSelfHighestFmaxRepa :: 
#   Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> 
#   Integer -> Integer ->
#   System -> Set.Set Variable -> HistoryRepa -> 
#   Maybe (System, DecompFud)

def parametersSystemsHistoryRepasDecomperMaxRollByMExcludedSelfHighestFmaxRepa(wmax,lmax,xmax,omax,bmax,mmax,umax,pmax,fmax,mult,seed,uu,vv,aa):
    repaRounding = 1e-6
    dom = relationsDomain
    def maxd(mm):
        if len(mm) > 0:
            return list(sset([(b,a) for (a,b) in mm]))[-1]
        return (0,sset())
    def tsgl(r):
        return sdict([(r,sdict())])
    uvars = systemsSetVar
    trim = histogramsTrim
    aall = histogramsList
    def red(aa,vv):
        return setVarsHistogramsReduce(vv,aa)
    def unit(ss):
        return setStatesHistogramUnit(sset([ss]))
    aahh = histogramsHistory
    hhhr = systemsHistoriesHistoryRepa
    def vars(hr):
        return sset(historyRepasVectorVar(hr))
    size = historyRepasSize
    rraa = systemsHistogramRepasHistogram
    hrhx = historyRepasRed
    def hrhrred(hr,vv):
        return setVarsHistoryRepasHistoryRepaReduced(vv,hr)
    def hrred(hr,vv):
        return setVarsHistoryRepasReduce(1,vv,hr)
    def reduce(uu,ww,hh):
        return rraa(uu,hrred(hh,ww))
    def select(uu,ss,hh):
        return historyRepasHistoryRepasHistoryRepaSelection_u(hhhr(uu,aahh(unit(ss))),hh)
    hrconcat = vectorHistoryRepasConcat_u
    hrshuffle = historyRepasShuffle_u
    ffqq = fudsSetTransform
    fder = fudsDerived
    tttr = systemsTransformsTransformRepa_u
    def apply(uu,ff,hh):
        return historyRepasListTransformRepasApply(hh,[tttr(uu,tt) for tt in ffqq(ff)])
    depends = fudsSetVarsDepends
    zzdf = treePairStateFudsDecompFud
    def zztrim(df):
        pp = []
        for ll in treesPaths(df):
            (_,ff) = ll[-1]
            if len(ff) == 0:
                pp.append(ll[:-1])
            else:
                pp.append(ll)
        return pathsTree(pp)
    def layerer(uu,xx,f):
        z = size(xx)
        xxp = hrhx(xx)
        xxrr = hrconcat([hrshuffle(xx,seed+i*z) for i in range(1,mult+1)])
        xxrrp = hrhx(xxrr)
        return parametersSystemsLayererMaxRollByMExcludedSelfHighestRepa(wmax,lmax,xmax,omax,bmax,mmax,umax,pmax,uu,vv,xx,xxp,xxrr,xxrrp,f)
    def decomp(uu,zz,qq,f):
        if len(zz) == 0:
            (uur,ffr,nnr) = layerer(uu,aa,f)
            if len(ffr) == 0 or len(nnr) == 0:
                return (uu, decompFudEmpty())
            (ar,kkr) = maxd(nnr)
            if ar <= repaRounding:
                return (uu, decompFudEmpty())
            ffr1 = depends(ffr,kkr)
            aar = apply(uur,ffr1,aa)
            aa1 = trim(reduce(uur,fder(ffr1),aar))
            zzr = tsgl((stateEmpty(),ffr1))
            qq[(stateEmpty(),ffr1)] = (aar,aa1)
            (ffr,nnr,kkr) = (None,None,None)
            return decomp(uur,zzr,qq,f+1)
        if fmax > 0 and f > fmax:
            return (uu,zzdf(zztrim(zz)))
        mm = []
        for (nn,yy) in treesPlaces(zz):
            (rr,ff) = nn[-1]
            if len(ff) > 0:
                (bb,bb1) = qq[(rr,ff)]
                tt = dom(treesRoots(yy))
                for (ss,a) in aall(red(bb1,fder(ff))):
                    if a > 0 and ss not in tt:
                        mm.append((a,(nn,ss,bb)))
        if len(mm) == 0:
            return (uu,zzdf(zztrim(zz)))
        mm.sort(key = lambda x: x[0])
        (_,(nn,ss,bb)) = mm[-1]
        cc = hrhrred(select(uu,ss,bb),vars(aa))
        (uuc,ffc,nnc) = layerer(uu,cc,f)
        (ac,kkc) = maxd(nnc)
        ffc1 = fudEmpty()
        if ac > repaRounding:
             ffc1 = depends(ffc,kkc)
        ccc = apply(uuc,ffc1,cc)
        cc1 = trim(reduce(uuc,fder(ffc1),ccc))
        qq[(ss,ffc1)] = (ccc,cc1)
        zzc = pathsTree(treesPaths(zz) + [nn+[(ss,ffc1)]])
        (mm,cc,ffc,nnc,kkc) = (None,None,None,None,None)
        return decomp(uuc,zzc,qq,f+1)
    if wmax < 0 or lmax < 0 or xmax < 0 or omax < 0 or bmax < 0 or mmax < 1 or umax < 0 or pmax < 0:
        return None
    if size(aa) == 0 or mult < 1:
        return None
    if not (vars(aa).issubset(uvars(uu)) and vv.issubset(vars(aa))):
        return None
    return decomp(uu,emptyTree(),sdict(),1)

# systemsDecompFudsHistoryRepasAlignmentContentShuffleSummation_u :: 
#   Integer -> Integer -> System -> DecompFud -> HistoryRepa -> (Double,Double)

def systemsDecompFudsHistoryRepasAlignmentContentShuffleSummation_u(mult,seed,uu,df,aa):
    vol = systemsSetVarsVolume_u
    vars = histogramsSetVar
    size = histogramsSize
    def resize(z,aa):
        if z > 0:
            return histogramsResize(z,aa)
        else:
            return histogramEmpty()
    algn = histogramsAlignment
    araa = systemsHistogramRepasHistogram
    def hrred(hr,vv):
        return setVarsHistoryRepasReduce(1,vv,hr)
    fder = fudsDerived
    apply = systemsDecompFudsHistoryRepasMultiplyWithShuffle
    qq = apply(mult,seed,uu,df,aa)
    a = 0.0
    ad = 0.0
    for ((_,ff),(hr,hrxx)) in qq.items():
        aa = araa(uu,hrred(hr,fder(ff)))
        u = vol(uu,vars(aa))
        m = len(vars(aa))
        bb = resize(size(aa),araa(uu,hrred(hrxx,fder(ff))))
        b = algn(aa) - algn(bb)
        a += b
        ad += b/(u**(1.0/m))
    return (a,ad)

# systemsDecompFudsHistoryRepasTreeAlignmentContentShuffleSummation_u :: 
#   Integer -> Integer -> System -> DecompFud -> HistoryRepa -> Map (State,Fud) (Int,(Double,Double))

def systemsDecompFudsHistoryRepasTreeAlignmentContentShuffleSummation_u(mult,seed,uu,df,aa):
    vol = systemsSetVarsVolume_u
    vars = histogramsSetVar
    size = histogramsSize
    def resize(z,aa):
        if z > 0:
            return histogramsResize(z,aa)
        else:
            return histogramEmpty()
    algn = histogramsAlignment
    araa = systemsHistogramRepasHistogram
    hrsize = historyRepasSize
    def hrred(hr,vv):
        return setVarsHistoryRepasReduce(1,vv,hr)
    fder = fudsDerived
    apply = systemsDecompFudsHistoryRepasMultiplyWithShuffle
    qq = apply(mult,seed,uu,df,aa)
    qq1 = sdict()
    for ((ss,ff),(hr,hrxx)) in qq.items():
        aa = araa(uu,hrred(hr,fder(ff)))
        u = vol(uu,vars(aa))
        m = len(vars(aa))
        bb = resize(size(aa),araa(uu,hrred(hrxx,fder(ff))))
        b = algn(aa) - algn(bb)
        qq1[(ss,ff)] = (hrsize(hr),(b,b/(u**(1.0/m))))
    return qq1

# parametersBuilderConditionalVarsRepa :: 
#   Integer -> Integer -> Integer -> Set.Set Variable -> HistoryRepa -> 
#   Maybe (Map.Map (Set.Set Variable) Double)

def parametersBuilderConditionalVarsRepa(kmax,omax,qmax,ll,aa):
    def sgl(x):
        return sset([x])
    def bot(amax,mm):
        return sdict([(x,y) for (y,x) in list(sset([(b,a) for (a,b) in mm.items()]))[:amax]])
    vars = historyRepasSetVariable
    def red(aa,vv):
        return setVarsHistoryRepasCountApproxs(vv,aa)
    if kmax < 0 or omax < 0 or qmax < 0:
        return None
    z = historyRepasSize(aa)
    def ent(xx):
        t = 0
        for i in xx:
            a = i/z
            t += a * log(a)
        return - t
    def buildc(qq,nn):
        pp = sset([kk|sgl(w) for (kk,e) in qq.items() if e > 0 for w in vvk-kk])
        mm = bot(omax,sdict([(jj,ent(red(aa,ll|jj))-ent(red(aa,jj))) for jj in pp if len(jj) <= kmax]))
        if len(mm) > 0:
            nn1 = nn.copy()
            nn1.update(mm)
            return buildc(mm,nn1)
        return nn
    vvk = vars(aa) - ll
    rr = bot(omax,sdict([(sgl(w),ent(red(aa,ll|sgl(w)))-ent(red(aa,sgl(w)))) for w in vvk]))
    return bot(qmax,buildc(rr,rr))

# parametersSystemsHistoryRepasDecomperConditionalFmaxRepa :: 
#   Integer -> Integer -> Integer -> System -> Set.Set Variable -> HistoryRepa -> 
#   Maybe (System, DecompFud)

def parametersSystemsHistoryRepasDecomperConditionalFmaxRepa(kmax,omax,fmax,uu,ll,aa):
    rounding = 1e-14
    dom = relationsDomain
    def tsgl(r):
        return sdict([(r,sdict())])
    uunion = pairSystemsUnion
    cart = systemsSetVarsSetStateCartesian_u
    llss = listsState
    sunion = pairStatesUnionLeft
    trim = histogramsTrim
    aall = histogramsList
    def red(aa,vv):
        return setVarsHistogramsReduce(vv,aa)
    unit = setStatesHistogramUnit
    mul = pairHistogramsMultiply
    ent = histogramsEntropy
    aahh = histogramsHistory
    hhhr = systemsHistoriesHistoryRepa
    def vars(hr):
        return sset(historyRepasVectorVar(hr))
    rraa = systemsHistogramRepasHistogram
    def hrhrred(hr,vv):
        return setVarsHistoryRepasHistoryRepaReduced(vv,hr)
    def hrred(hr,vv):
        return setVarsHistoryRepasReduce(1,vv,hr)
    def reduce(uu,ww,hh):
        return rraa(uu,hrred(hh,ww))
    def select(uu,ss,hh):
        return historyRepasHistoryRepasHistoryRepaSelection_u(hhhr(uu,aahh(unit(sset([ss])))),hh)
    trans = histogramsSetVarsTransform_u
    tttr = systemsTransformsTransformRepa_u
    fsys = fudsSystemImplied
    qqff = setTransformsFud_u
    ffqq = fudsSetTransform
    fder = fudsDerived
    def apply(uu,ff,hh):
        return historyRepasListTransformRepasApply(hh,[tttr(uu,tt) for tt in ffqq(ff)])
    zzdf = treePairStateFudsDecompFud
    def vvff(uu,vv,f):
        v = VarPair((VarPair((VarInt(f),VarInt(1))),VarInt(1)))
        qq = sset([sunion(ss,llss([(v,ValInt(i+1))])) for (i,ss) in enumerate(cart(uu,vv))])
        return qqff(sset([trans(unit(qq),sset([v]))]))
    def lenter(ll,aa):
        return list(parametersBuilderConditionalVarsRepa(kmax,omax,1,ll,aa).items())
    vv = vars(aa)
    def decomp(uu,zz,qq,f):
        if len(zz) == 0:
            nnr = lenter(ll,aa)
            if len(nnr) == 0:
                return (uu, decompFudEmpty())
            [(kkr,_)] = nnr
            ffr = vvff(uu,kkr,f)
            uur = uunion(uu,fsys(ffr))
            aar = apply(uur,ffr,aa)
            aa1 = trim(reduce(uur,fder(ffr)|ll,aar))
            zzr = tsgl((stateEmpty(),ffr))
            qq[(stateEmpty(),ffr)] = (aar,aa1)
            (ffr,nnr,kkr) = (None,None,None)
            return decomp(uur,zzr,qq,f+1)
        if fmax > 0 and f > fmax:
            return (uu,zzdf(zz))
        mm = []
        for (nn,yy) in treesPlaces(zz):
            (rr,ff) = nn[-1]
            (bb,bb1) = qq[(rr,ff)]
            tt = dom(treesRoots(yy))
            for (ss,a) in aall(red(bb1,fder(ff))):
                if a > 0 and ss not in tt:
                    e = a * ent(red(mul(bb1,unit(sset([ss]))),ll))
                    if e > rounding:
                        mm.append((e,(nn,ss,bb)))
        if len(mm) == 0:
            return (uu,zzdf(zz))
        mm.sort(key = lambda x: x[0])
        (_,(nn,ss,bb)) = mm[-1]
        cc = hrhrred(select(uu,ss,bb),vv)
        nnc = lenter(ll,cc)
        [(kkc,_)] = nnc
        ffc = vvff(uu,kkc,f)
        uuc = uunion(uu,fsys(ffc))
        ccc = apply(uuc,ffc,cc)
        cc1 = trim(reduce(uuc,fder(ffc)|ll,ccc))
        qq[(ss,ffc)] = (ccc,cc1)
        zzc = pathsTree(treesPaths(zz) + [nn+[(ss,ffc)]])
        (mm,cc,ffc,nnc,kkc) = (None,None,None,None,None)
        return decomp(uuc,zzc,qq,f+1)
    if kmax < 0 or omax < 0:
        return None
    return decomp(uu,emptyTree(),sdict(),1)