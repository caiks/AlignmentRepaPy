from AlignmentApprox import *
import numpy as np
import scipy.special as sp

AlignmentForeignPy_ok = True
try:
    from AlignmentForeignPy import *
except ImportError:
    AlignmentForeignPy_ok = False

# data HistogramRepa = HistogramRepa {
#   histogramRepasVectorVar :: !(V.Vector Variable),
#   histogramRepasMapVarInt :: Map.Map Variable Int,
#   histogramRepasArray :: !(Array U VShape Double)}

# histogramRepaEmpty :: HistogramRepa

def histogramRepaEmpty():
    return ([],sdict(), np.empty(0,dtype=np.dtype(np.float64)))

# histogramRepasVectorVar :: HistogramRepa -> V.Vector Variable

def histogramRepasVectorVar(ar):
    (vv,_,_) = ar
    return vv

# histogramRepasMapVarInt :: HistogramRepa -> Map.Map Variable Int

def histogramRepasMapVarInt(ar):
    (_,mvv,_) = ar
    return mvv

# histogramRepasArray :: HistogramRepa -> Array U VShape Double

def histogramRepasArray(ar):
    (_,_,rr) = ar
    return rr

# arraysHistogramRepa :: Array U VShape Double -> HistogramRepa

def arraysHistogramRepa(rr):
    n = rr.ndim
    vv = [VarIndex(i) for i in range(0,n)]
    mvv = sdict([(v,i) for (i,v) in enumerate(vv)])
    return (vv,mvv,rr)

# systemsHistogramsHistogramRepa :: System -> Histogram -> Maybe HistogramRepa

def systemsHistogramsHistogramRepa(uu,aa):
    uvals = systemsVarsSetValue
    uvars = systemsSetVar
    sat = statesVarsValue
    vars = histogramsSetVar
    aall = histogramsList
    if len(aa) > 0 and vars(aa).issubset(uvars(uu)):
        vv = list(vars(aa))
        mvv = sdict([(v,i) for (i,v) in enumerate(vv)])
        mm = sdict([(v,sdict([(w,i) for (i,w) in enumerate(uvals(uu,v))])) for v in vv])
        sh = tuple([len(mm[v]) for v in vv])
        rr = np.zeros(sh,dtype=np.dtype(np.float64))
        for (ss,c) in aall(aa):
            rr[tuple([mm[v][sat(ss,v)] for v in vv])] = float(c)
        return (vv,mvv,rr)
    return None

# systemsHistogramRepasHistogram :: System -> HistogramRepa -> Maybe Histogram

def systemsHistogramRepasHistogram(uu,ar):
    uvals = systemsVarsSetValue
    uvars = systemsSetVar
    llss = listsState
    llaa = listsHistogram
    (vv,_,rr) = ar
    mm = sdict([(v,list(uvals(uu,v))) for v in vv])
    sh = tuple([len(mm[v]) for v in vv])
    sh1 = rr.shape
    if sset(vv).issubset(uvars(uu)) and sh == sh1:
        return llaa([(llss([(vv[j],mm[vv[j]][k]) for (j,k) in enumerate(ii)]), ratio(c)) for (ii,c) in np.ndenumerate(rr)])
    return None

# setVarsHistogramRepasReduce :: Set.Set Variable -> HistogramRepa -> HistogramRepa 

def setVarsHistogramRepasReduce(kk,ar):
    (vvv,mvv,rr) = ar
    vv = sset(vvv)
    if len(vv - kk) == 0:
        return ar
    vkk = list(kk & vv)
    if len(vkk) == 0:
        return ([],sdict(),np.sum(rr))
    mkk = sdict([(v,i) for (i,v) in enumerate(vkk)])
    pkk = tuple([mvv[v] for v in vv-vkk])
    return (vkk,mkk,np.sum(rr,axis=pkk))

# setSetVarsHistogramRepasPartitionIndependent_u :: Set.Set (Set.Set Variable) -> Double -> HistogramRepa -> HistogramRepa 

def setSetVarsHistogramRepasPartitionIndependent_u(pp,z,ar):
    def prod(ll):
        p = 1
        for x in ll:
            p *= x
        return p
    (vvv,mvv,rr) = ar
    n = len(vvv)
    m = len(pp)
    f = 1.0 / z
    svv = rr.shape
    ppp = [[mvv[v] for v in cc] for cc in pp]
    sxx = tuple([prod([svv[p] for p in pcc]) for pcc in ppp])
    pww = [p for pcc in ppp for p in pcc]
    rr1 = np.reshape(np.moveaxis(rr,pww,list(range(n))),sxx)
    rr2 = np.sum(rr1,axis=tuple(range(1,m)))
    for i in range(1,m):
        rr2 = np.multiply.outer(rr2,np.sum(rr1,axis=tuple([j for j in range(m) if j != i]))*f)
    return arraysHistogramRepa(rr2)

# sumsHistogramRepasRollMapPairsHistogramRepasSum_u :: Double -> HistogramRepa -> 
#   (UV.Vector Int, UV.Vector Int) -> HistogramRepa -> UV.Vector Double

def sumsHistogramRepasRollMapPairsHistogramRepasSum_u(a,aav,xx,rrv):
    (_,_,aa) = aav
    (_,_,rr) = rrv
    (rs,rt) = xx
    return rr + a - aa[rs] - aa[rt]


# data HistogramRepaVec = HistogramRepaVec {
#   histogramRepaVecsVectorVar :: !(V.Vector Variable),
#   histogramRepaVecsMapVarInt :: Map.Map Variable Int,
#   histogramRepaVecsSize :: !Double,
#   histogramRepaVecsShape :: !VShape,
#   histogramRepaVecsArray :: !(V.Vector (UV.Vector Double))}

# histogramRepaVecEmpty :: HistogramRepaVec

def histogramRepaVecEmpty():
    return ([],sdict(),0,tuple(),[])

# histogramRepaVecsVectorVar :: HistogramRepaVec -> V.Vector Variable

def histogramRepaVecsVectorVar(hr):
    return hr[0]

# histogramRepaVecsMapVarInt :: HistogramRepaVec -> Map.Map Variable Int

def histogramRepaVecsMapVarInt(hr):
    return hr[1]

# histogramRepaVecsSize :: HistogramRepaVec -> V.Double

def histogramRepaVecsSize(hr):
    return hr[2]

# histogramRepaVecsShape :: HistogramRepaVec -> VShape

def histogramRepaVecsShape(hr):
    return hr[3]

# histogramRepaVecsArray :: HistogramRepaVec -> V.Vector (UV.Vector Double)

def histogramRepaVecsArray(hr):
    return hr[4]

# vectorHistogramRepasHistogramRepaVec_u :: Double -> V.Vector HistogramRepa -> HistogramRepaVec

def vectorHistogramRepasHistogramRepaVec_u(z,vrr):
    rraa = histogramRepasArray
    (vvv,mvv,rr) = vrr[0]
    svv = rr.shape
    return (vvv,mvv,z,svv,list(map(rraa,vrr)))

# histogramRepaVecsSum :: HistogramRepaVec -> UV.Vector Double

def histogramRepaVecsSum(vhr):
    (_,_,_,_,vaa) = vhr
    return list(map(np.sum,vaa))

# histogramRepaVecsFaclnsRepaVecs :: HistogramRepaVec -> HistogramRepaVec

def histogramRepaVecsFaclnsRepaVecs(vhr):
    (vvv,mvv,_,svv,vaa) = vhr
    return (vvv,mvv,1,svv,[sp.gammaln(aa + 1) for aa in vaa])

# setSetVarsHistogramRepaVecsPartitionVec_u :: Set.Set (Set.Set Variable) -> HistogramRepaVec -> HistogramRepaVec 

def setSetVarsHistogramRepaVecsPartitionVec_u(pp,vhr):
    def prod(ll):
        p = 1
        for x in ll:
            p *= x
        return p
    (vvv,mvv,z,svv,vaa) = vhr
    n = len(vvv)
    m = len(pp)
    ppp = [[mvv[v] for v in cc] for cc in pp]
    svv1 = tuple([prod([svv[p] for p in pcc]) for pcc in ppp])
    pww0 = list(range(n))
    pww = [p for pcc in ppp for p in pcc]
    vaa1 = []
    for rr in vaa:
        if pww != pww0:
            rr0 = np.moveaxis(rr,pww,pww0)
        else:
            rr0 = rr
        rr1 = np.reshape(rr0,svv1)
        vaa1.append(rr1)
    vaa1 = [np.reshape(np.moveaxis(rr,pww,list(range(n))),svv1) for rr in vaa]
    vvv1 = [VarIndex(i) for i in range(m)]
    mvv1 = sdict([(v,i) for (i,v) in enumerate(vvv1)])
    return (vvv1,mvv1,z,svv1,vaa1)

# setSetVarsHistogramRepaVecsPartitionVec_u_2 :: Set.Set (Set.Set Variable) -> HistogramRepaVec -> HistogramRepaVec 

def setSetVarsHistogramRepaVecsPartitionVec_u_2(pp,vhr):
    def prod(ll):
        p = 1
        for x in ll:
            p *= x
        return p
    (vvv,mvv,z,svv,vaa) = vhr
    n = len(vvv)
    m = len(pp)
    ppp = [[mvv[v] for v in cc] for cc in pp]
    svv1 = tuple([prod([svv[p] for p in pcc]) for pcc in ppp])
    pww = [p for pcc in ppp for p in pcc]
    vaa1 = [np.reshape(np.moveaxis(rr,pww,list(range(n))),svv1) for rr in vaa]
    vvv1 = [VarIndex(i) for i in range(m)]
    mvv1 = sdict([(v,i) for (i,v) in enumerate(vvv1)])
    return (vvv1,mvv1,z,svv1,vaa1)

# setSetVarsHistogramRepaVecsPartitionIndependentVec_u :: Set.Set (Set.Set Variable) -> HistogramRepaVec -> HistogramRepaVec 

def setSetVarsHistogramRepaVecsPartitionIndependentVec_u(pp,vhr):
    def prod(ll):
        p = 1
        for x in ll:
            p *= x
        return p
    (vvv,mvv,z,svv,vaa) = vhr
    n = len(vvv)
    m = len(pp)
    f = 1.0 / z
    ppp = [[mvv[v] for v in cc] for cc in pp]
    svv1 = tuple([prod([svv[p] for p in pcc]) for pcc in ppp])
    pww0 = list(range(n))
    pww = [p for pcc in ppp for p in pcc]
    vaa1 = []
    for rr in vaa:
        if pww != pww0:
            rr0 = np.moveaxis(rr,pww,pww0)
        else:
            rr0 = rr
        rr1 = np.reshape(rr0,svv1)
        rr2 = np.sum(rr1,axis=tuple(range(1,m)))
        for i in range(1,m):
            rr2 = np.multiply.outer(rr2,np.sum(rr1,axis=tuple([j for j in range(m) if j != i]))*f)
        vaa1.append(rr2)
    vvv1 = [VarIndex(i) for i in range(m)]
    mvv1 = sdict([(v,i) for (i,v) in enumerate(vvv1)])
    return (vvv1,mvv1,z,svv1,vaa1)

# setSetVarsHistogramRepaVecsPartitionIndependentVec_u_1 :: Set.Set (Set.Set Variable) -> HistogramRepaVec -> HistogramRepaVec 

def setSetVarsHistogramRepaVecsPartitionIndependentVec_u_1(pp,vhr):
    def prod(ll):
        p = 1
        for x in ll:
            p *= x
        return p
    (vvv,mvv,z,svv,vaa) = vhr
    n = len(vvv)
    m = len(pp)
    f = 1.0 / z
    ppp = [[mvv[v] for v in cc] for cc in pp]
    svv1 = tuple([prod([svv[p] for p in pcc]) for pcc in ppp])
    pww = [p for pcc in ppp for p in pcc]
    vaa1 = []
    for rr in vaa:
        rr1 = np.reshape(np.moveaxis(rr,pww,list(range(n))),svv1)
        rr2 = np.sum(rr1,axis=tuple(range(1,m)))
        for i in range(1,m):
            rr2 = np.multiply.outer(rr2,np.sum(rr1,axis=tuple([j for j in range(m) if j != i]))*f)
        vaa1.append(rr2)
    vvv1 = [VarIndex(i) for i in range(m)]
    mvv1 = sdict([(v,i) for (i,v) in enumerate(vvv1)])
    return (vvv1,mvv1,z,svv1,vaa1)

# varsHistogramRepa4VecsReduceSingle_u :: Int -> HistogramRepaVec -> HistogramRepa

def varsHistogramRepa4VecsReduceSingle_u(u,vhr):
    (vvv,_,_,_,vaa) = vhr
    n = len(vvv)
    x = vvv[u]
    vyy = [x]
    myy = sdict([(x,0)])
    ax = tuple([j for j in range(n) if j != u])
    a1 = np.sum(vaa[0],axis=ax)
    a2 = np.sum(vaa[1],axis=ax)
    b1 = np.sum(vaa[2],axis=ax)
    b2 = np.sum(vaa[3],axis=ax)
    return (vyy,myy,a1-a2-b1+b2)

# varsHistogramRepaVecsRollVec_u :: Int -> HistogramRepaVec -> (HistogramRepaVec, (UV.Vector Int, UV.Vector Int))

def varsHistogramRepaVecsRollVec_u(u,vhr):
    (vvv,mvv,_,svv,vaa) = vhr
    p = len(vaa)
    d = svv[u]
    rr = list(zip(*[(s,t) for s in range(1,d) for t in range(d-1) if s > t]))
    rs = list(rr[0])
    rt = list(rr[1])
    r = len(rs)
    syy0 = list(svv)
    syy0[u] = r
    syy = tuple(syy0)
    vbb = []
    for i in range(p):
        aa = np.swapaxes(vaa[i],0,u)
        vbb.append(np.swapaxes(aa[rs]+aa[rt],0,u))
    return ((vvv,mvv,1,syy,vbb),(rs,rt))

# varsSourcesTargetsRollsHistogramRepaVecsHistogramRepaVecRollsCopyVec_u :: 
#   Int -> Int -> Int -> Int -> HistogramRepaVec -> HistogramRepaVec -> HistogramRepaVec

def varsSourcesTargetsRollsHistogramRepaVecsHistogramRepaVecRollsCopyVec_u(u,s,t,q,rrv,ssv):
    (vvv,mvv,z,svv,vaa) = rrv
    (_,_,_,syy,vyy) = ssv
    p = len(vaa)
    d = svv[u]
    sbb0 = list(svv)
    sbb0[u] = d-1
    sbb = tuple(sbb0)
    vbb = []
    for i in range(p):
        aa = np.swapaxes(vaa[i],0,u)
        yy = np.swapaxes(vyy[i],0,u)
        vbb.append(np.swapaxes(np.concatenate((aa[:t],yy[[q]],aa[t+1:s],aa[s+1:])),0,u))
    return (vvv,mvv,z,sbb,vbb)


# data HistogramRepaRed = HistogramRepaRed {
#   histogramRepaRedsVectorVar :: !(V.Vector Variable),
#   histogramRepaRedsMapVarInt :: Map.Map Variable Int,
#   histogramRepaRedsShape :: !VShape,
#   histogramRepaRedsVectorArray :: !(V.Vector (UV.Vector Double))}
#                deriving (Eq, Read, Show)

# histogramRepaRedsVectorVar :: HistogramRepaRed -> V.Vector Variable

def histogramRepaRedsVectorVar(arr):
    (vv,_,_,_) = arr
    return vv

# histogramRepaRedsMapVarInt :: HistogramRepaRed -> Map.Map Variable Int

def histogramRepaRedsMapVarInt(arr):
    (_,mvv,_,_) = arr
    return mvv

# histogramRepaRedsShape :: HistogramRepaRed -> VShape

def histogramRepaRedsShape(arr):
    (_,_,sh,_) = arr
    return sh

# histogramRepaRedsVectorArray :: HistogramRepaRed -> Array U VShape Double

def histogramRepaRedsVectorArray(arr):
    (_,_,_,lrr) = arr
    return lrr

# histogramRepasRed_u :: Double -> HistogramRepa -> HistogramRepaRed 

def histogramRepasRed_u(z,ar):
    (vvv,mvv,rr) = ar
    f = 1.0 / z
    n = len(vvv)
    lrr = []
    for v in vvv:
        i = mvv[v]
        lrr.append(np.sum(rr,axis=tuple([j for j in range(n) if j != i])))
    if f != 1.0:
        lrr = [ss * f for ss in lrr]
    return (vvv,mvv,rr.shape,lrr)

# histogramRepaRedsIndependent :: Double -> HistogramRepaRed -> HistogramRepa 

def histogramRepaRedsIndependent(z,arr):
    (vvv,mvv,svv,lrr) = arr
    if len(vvv) == 0:
        return ([],sdict(),np.array(z,dtype=np.dtype(np.float64)))
    rr = lrr[0]
    for i in range(1,len(vvv)):
        rr = np.multiply.outer(rr,lrr[i])
    return (vvv,mvv,rr*z)

# setVarsHistogramRepaRedsRed :: Set.Set Variable -> HistogramRepaRed -> HistogramRepaRed 

def setVarsHistogramRepaRedsRed(kk,arr):
    (vvv,mvv,svv,lrr) = arr
    vv = sset(vvv)
    vkk = list(kk & vv)
    if len(vkk) == 0:
        return ([],sdict(),tuple(),[])
    mkk = sdict([(v,i) for (i,v) in enumerate(vkk)])
    pkk = [mvv[v] for v in vkk]
    skk = tuple([svv[i] for i in pkk])
    lkk = [lrr[i] for i in pkk]
    return (vkk,mkk,skk,lkk)

# varsHistogramRepaRedsSingle_u :: Int -> HistogramRepaRed -> HistogramRepa

def varsHistogramRepaRedsSingle_u(u,xxv):
    (vvv,_,_,vxx) = xxv
    x = vvv[u]
    vyy = [x]
    myy = sdict([(x,0)])
    bb = vxx[u]
    return (vyy,myy,bb)

# histogramRepa4VecsRed_u :: HistogramRepaVec -> HistogramRepaRed

def histogramRepa4VecsRed_u(vhr):
    (vvv,mvv,_,svv,vaa) = vhr
    n = len(vvv)
    vbb = []
    for u in range(n):
        ax = tuple([j for j in range(n) if j != u])
        a1 = np.sum(vaa[0],axis=ax)
        a2 = np.sum(vaa[1],axis=ax)
        b1 = np.sum(vaa[2],axis=ax)
        b2 = np.sum(vaa[3],axis=ax)
        vbb.append(a1-a2-b1+b2)
    return (vvv,mvv,svv,vbb)


# data HistoryRepa = HistoryRepa {
#   historyRepasVectorVar :: !(V.Vector Variable),
#   historyRepasMapVarInt :: Map.Map Variable Int,
#   historyRepasShape :: !VShape,
#   historyRepasArray :: !(Array U DIM2 Int)}

# historyRepaEmpty :: HistoryRepa

def historyRepaEmpty():
    return ([],sdict(),tuple(),np.empty((0,0),dtype=np.dtype(np.int32)))

# historyRepasVectorVar :: HistoryRepa -> V.Vector Variable

def historyRepasVectorVar(hr):
    (vv,_,_,_) = hr
    return vv

# historyRepasMapVarInt :: HistoryRepa -> Map.Map Variable Int

def historyRepasMapVarInt(hr):
    (_,mvv,_,_) = hr
    return mvv

# historyRepasShape :: HistoryRepa -> VShape

def historyRepasShape(hr):
    (_,_,sh,_) = hr
    return sh

# historyRepasArray :: HistoryRepa -> Array U DIM2 Int

def historyRepasArray(hr):
    (_,_,sh,rr) = hr
    return rr

# systemsHistoriesHistoryRepa :: System -> History -> Maybe HistoryRepa

def systemsHistoriesHistoryRepa(uu,hh):
    uvals = systemsVarsSetValue
    uvars = systemsSetVar
    vars = historiesSetVar
    ssll = statesList
    hhll = historiesList
    if len(hh) > 0 and vars(hh).issubset(uvars(uu)):
        vv = list(vars(hh))
        mvv = sdict([(v,i) for (i,v) in enumerate(vv)])
        mm = sdict([(v,sdict([(w,i) for (i,w) in enumerate(uvals(uu,v))])) for v in vv])
        sh = tuple([len(mm[v]) for v in vv])
        sh1 =(len(hh),len(vv))
        nn = [mm[v][u] for (_,ss) in hhll(hh) for (v,u) in ssll(ss)]
        rr = np.transpose(np.reshape(np.array(nn,dtype=np.dtype(np.int32)),sh1))
        return (vv,mvv,sh,rr)
    return None

# systemsHistoryRepasHistory_u :: System -> HistoryRepa -> Maybe History

def systemsHistoryRepasHistory_u(uu,hr):
    uvals = systemsVarsSetValue
    uvars = systemsSetVar
    llss = listsState
    llhh = listsHistory
    (vv,_,_,rr) = hr
    mm = sdict([(v,list(uvals(uu,v))) for v in vv])
    (n,z) = rr.shape
    if sset(vv).issubset(uvars(uu)):
        return llhh([(IdInt(i+1),llss([(vv[j],mm[vv[j]][rr[j,i]]) for j in range(0,n)])) for i in range(0,z)])
    return None

# vectorHistoryRepasConcat_u :: V.Vector HistoryRepa -> HistoryRepa 

def vectorHistoryRepasConcat_u(ll):
    (vaa,maa,saa,_) = ll[0]
    rbb = np.transpose(np.concatenate([np.transpose(rr) for (_,_,_,rr) in ll]))
    return (vaa,maa,saa,rbb)

# eventsHistoryRepasHistoryRepaSelection :: [Int] -> HistoryRepa -> HistoryRepa 

def eventsHistoryRepasHistoryRepaSelection(ll,hr):
    (vvv,mvv,svv,rr) = hr
    (n,z) = rr.shape
    ll1 = [i for i in ll if i >= 0 and i < z]
    rr1 = np.transpose(np.transpose(rr)[ll1])
    return (vvv,mvv,svv,rr1)

# historyRepasHistoryRepasHistoryRepaSelection_u :: HistoryRepa -> HistoryRepa -> HistoryRepa 

def historyRepasHistoryRepasHistoryRepaSelection_u(ss,hh):
    if ss == historyRepaEmpty() or hh == historyRepaEmpty():
        return historyRepaEmpty()
    (vss,_,_,rss) = ss
    (vhh,mhh,shh,rhh) = hh
    if len(vss) == 0:
        return hh
    (_,y) = rss.shape
    (_,z) = rhh.shape
    pss = [mhh[v] for v in vss]
    rhh0 = np.transpose(rhh[pss])
    rss0 = np.transpose(rss)
    ll = []
    for i in range(z):
        aa = rhh0[i]
        for j in range(y):
            if np.all(rss0[j]==aa):
                ll.append(i)
                break
    rhh1 = np.transpose(np.transpose(rhh)[ll])
    return (vhh,mhh,shh,rhh1)



# setVarsHistoryRepasHistoryRepaReduced :: Set.Set Variable -> HistoryRepa -> HistoryRepa 

def setVarsHistoryRepasHistoryRepaReduced(kk,hr):
    (vvv,mvv,svv,rr) = hr
    vv = sset(vvv)
    vkk = list(kk & vv)
    mkk = sdict([(v,i) for (i,v) in enumerate(vkk)])
    pkk = [mvv[v] for v in vkk]
    skk = tuple([svv[i] for i in pkk])
    rr1 = rr[pkk]
    return (vkk,mkk,skk,rr1)

# void listVarsArrayHistoriesReduce_u(double f, int n, int* ppkk, int* pskk, int z, int* prr, double* pmv)

# setVarsHistoryRepasReduce :: Double -> Set.Set Variable -> HistoryRepa -> HistogramRepa 

def setVarsHistoryRepasReduce(f,kk,hr):
    (vvv,mvv,svv,rr) = hr
    vv = sset(vvv)
    (_,z) = rr.shape
    vkk = list(kk & vv)
    m = len(vkk)
    if m == 0:
        return ([],sdict(),np.array(z,dtype=np.dtype(np.float64)))
    mkk = sdict([(v,i) for (i,v) in enumerate(vkk)])
    pkk = [mvv[v] for v in vkk]
    skk = tuple([svv[i] for i in pkk])
    rr1 = np.zeros(skk,dtype=np.dtype(np.float64))
    listVarsArrayHistoriesReduce_u_py(f,m,np.array(pkk,dtype=np.dtype(np.int32)),np.array(skk,dtype=np.dtype(np.int32)),z,rr,rr1)
    return (vkk,mkk,rr1)

# setVarsHistoryRepasReduce_1 :: Double -> Set.Set Variable -> HistoryRepa -> HistogramRepa 

def setVarsHistoryRepasReduce_1(f,kk,hr):
    (vvv,mvv,svv,rr) = hr
    vv = sset(vvv)
    (_,z) = rr.shape
    vkk = list(kk & vv)
    if len(vkk) == 0:
        return ([],sdict(),np.array(z,dtype=np.dtype(np.float64)))
    mkk = sdict([(v,i) for (i,v) in enumerate(vkk)])
    pkk = [mvv[v] for v in vkk]
    skk = tuple([svv[i] for i in pkk])
    rr1 = np.zeros(skk,dtype=np.dtype(np.float64))
    rrt = np.transpose(rr)
    for j in range(0,z):
        rr1[tuple(rrt[j][pkk])] += f
    return (vkk,mkk,rr1)

# setVarsHistoryRepasReduce_2 :: Double -> Set.Set Variable -> HistoryRepa -> HistogramRepa 

def setVarsHistoryRepasReduce_2(f,kk,hr):
    (vvv,mvv,svv,rr) = hr
    vv = sset(vvv)
    (_,z) = rr.shape
    vkk = list(kk & vv)
    if len(vkk) == 0:
        return ([],sdict(),np.array(z,dtype=np.dtype(np.float64)))
    mkk = sdict([(v,i) for (i,v) in enumerate(vkk)])
    pkk = [mvv[v] for v in vkk]
    skk = tuple([svv[i] for i in pkk])
    rr1 = np.zeros(skk,dtype=np.dtype(np.float64))
    rrkr = np.transpose(rr[pkk])
    for j in range(0,z):
        rr1[tuple(rrkr[j])] += f
    return (vkk,mkk,rr1)

# setVarsHistoryRepasReduce_3 :: Double -> Set.Set Variable -> HistoryRepa -> HistogramRepa 

def setVarsHistoryRepasReduce_3(f,kk,hr):
    (vvv,mvv,svv,rr) = hr
    vv = sset(vvv)
    (_,z) = rr.shape
    vkk = list(kk & vv)
    m = len(vkk)
    if m == 0:
        return ([],sdict(),np.array(z,dtype=np.dtype(np.float64)))
    mkk = sdict([(v,i) for (i,v) in enumerate(vkk)])
    pkk = [mvv[v] for v in vkk]
    skk = tuple([svv[i] for i in pkk])
    rr1 = np.zeros(tuple(list(skk)+[z]),dtype=np.dtype(np.float64))
    rrkr = np.concatenate([rr[pkk],np.reshape(np.arange(z),(1,z))])
    rr1[tuple([rrkr[i] for i in range(m+1)])] = f
    return (vkk,mkk,np.sum(rr1,axis=-1))

if not AlignmentForeignPy_ok:
    setVarsHistoryRepasReduce = setVarsHistoryRepasReduce_3


# historyRepasRed :: HistoryRepa -> HistogramRepaRed 

def historyRepasRed(hr):
    (vvv,mvv,svv,rr) = hr
    (n,z) = rr.shape
    f = 1.0 / z
    lrr = [np.zeros(i,dtype=np.dtype(np.float64)) for i in svv]
    for i in range(0,n):
        for w in np.nditer(rr[i]):
            lrr[i][w] += f
    return (vvv,mvv,svv,lrr)

# setVarsHistoryRepasRed :: Set.Set Variable -> HistoryRepa -> HistogramRepaRed 

def setVarsHistoryRepasRed(kk,hr):
    (vvv,mvv,svv,rr) = hr
    vv = sset(vvv)
    vkk = list(kk & vv)
    if len(vkk) == 0:
        return ([],sdict(),tuple(),[])
    (_,z) = rr.shape
    mkk = sdict([(v,i) for (i,v) in enumerate(vkk)])
    pkk = [mvv[v] for v in vkk]
    skk = tuple([svv[i] for i in pkk])
    f = 1.0 / z
    lrr = [np.zeros(i,dtype=np.dtype(np.float64)) for i in skk]
    for (i,j) in enumerate(pkk):
        for w in np.nditer(rr[j]):
            lrr[i][w] += f
    return (vkk,mkk,skk,lrr)

# historyRepasSize :: HistoryRepa -> Int 

def historyRepasSize(hr):
    (_,_,_,rr) = hr
    (_,z) = rr.shape
    return z

# historyRepasSetVariable :: HistoryRepa -> Set.Set Variable 

def historyRepasSetVariable(hr):
    (vvv,_,_,_) = hr
    return sset(vvv)

# setVarsHistoryRepasCountApproxs :: Set.Set Variable -> HistoryRepa -> [Int] 

def setVarsHistoryRepasCountApproxs(kk,hr):
    (vvv,mvv,_,rr) = hr
    vv = sset(vvv)
    vkk = list(kk & vv)
    m = len(vkk)
    (_,z) = rr.shape
    if m == 0:
        return [z]
    pkk = [mvv[v] for v in vkk]
    rrk = rr[pkk]
    rrk1 = np.copy(rrk[0])
    for j in range(1,m):
        rrk1 *= 23 
        rrk1 += rrk[j]
    ll = dict()
    for a in rrk1.flat:
        ll[a] = ll.get(a,0) + 1
    return list(ll.values())

# setVarsHistoryRepasCountApproxs_1 :: Set.Set Variable -> HistoryRepa -> [Int] 

def setVarsHistoryRepasCountApproxs_1(kk,hr):
    (vvv,mvv,_,rr) = hr
    vv = sset(vvv)
    vkk = list(kk & vv)
    (_,z) = rr.shape
    if len(vkk) == 0:
        return [z]
    pkk = [mvv[v] for v in vkk]
    rrt = np.transpose(rr)
    ll = dict()
    for j in range(0,z):
        a = 0
        for k in np.nditer(rrt[j][pkk]):
            a = a * 23 + k
        ll[a] = ll.get(a,0) + 1
    return list(ll.values())



# data TransformRepa = TransformRepa {
#   transformRepasVectorVar :: !(V.Vector Variable),
#   transformRepasMapVarInt :: Map.Map Variable Int,
#   transformRepasVarDerived :: !Variable,
#   transformRepasValency :: !Int,
#   transformRepasArray :: !(Array U VShape Int)}


# transformRepasVectorVar :: TransformRepa -> V.Vector Variable

def transformRepasVectorVar(tr):
    (vv,_,_,_,_) = tr
    return vv

# transformRepasMapVarInt :: TransformRepa -> Map.Map Variable Int

def transformRepasMapVarInt(tr):
    (_,mvv,_,_,_) = tr
    return mvv

# transformRepasVarDerived :: TransformRepa -> Variable

def transformRepasVarDerived(tr):
    (_,_,w,_,_) = tr
    return w

# transformRepasValency :: TransformRepa -> Int

def transformRepasValency(tr):
    (_,_,_,d,_) = tr
    return d

# transformRepasArray :: TransformRepa -> Array U VShape Int

def transformRepasArray(tr):
    (_,_,_,_,rr) = tr
    return rr

# systemsTransformsTransformRepa_u :: System -> Transform -> TransformRepa

def systemsTransformsTransformRepa_u(uu,tt):
    uvals = systemsVarsSetValue
    sat = statesVarsValue
    und = transformsUnderlying
    der = transformsDerived
    aall = histogramsList
    ttaa = transformsHistogram
    vv = list(und(tt))
    w = der(tt)[0]
    mvv = sdict([(v,i) for (i,v) in enumerate(vv)])
    mm = sdict([(v,sdict([(w,i) for (i,w) in enumerate(uvals(uu,v))])) for v in vv + [w]])
    sh = tuple([len(mm[v]) for v in vv])
    d = len(mm[w])
    rr = np.empty(sh,dtype=np.dtype(np.int32))
    for (ss,_) in aall(ttaa(tt)):
        rr[tuple([mm[v][sat(ss,v)] for v in vv])] = mm[w][sat(ss,w)]
    return (vv,mvv,w,d,rr)

# historyRepasTransformRepasApply_u :: HistoryRepa -> TransformRepa -> HistoryRepa 

def historyRepasTransformRepasApply_u(hr,tr):
    (_,maa,_,raa) = hr
    (vtt,_,w,d,rtt) = tr
    (_,z) = raa.shape
    pkk = [maa[v] for v in vtt]
    rbb = np.empty(z,dtype=np.dtype(np.int32))
    rat = np.transpose(raa)
    for j in range(0,z):
        rbb[j] = rtt[tuple(rat[j][pkk])]
    return ([w],sdict([(w,0)]),(d,),np.reshape(rbb,(1,z)))

# historyRepasListTransformRepasApply_u :: HistoryRepa -> V.Vector TransformRepa -> HistoryRepa 

def historyRepasListTransformRepasApply_u(hr,fr):
    (vaa,maa,saa,raa) = hr
    (n,z) = raa.shape
    m = len(fr)
    rbt = np.transpose(np.append(raa,np.empty((m,z),dtype=np.dtype(np.int32)),axis=0))
    vbb = vaa
    sbb = list(saa)
    mbb = maa.copy()
    for (q,tr) in enumerate(fr):
        (vtt,_,w,d,rtt) = tr
        vbb = vbb + [w]
        sbb = sbb + [d]
        i = n + q
        mbb[w] = i
        pkk = [mbb[v] for v in vtt]
        for j in range(0,z):
            rbt[j,i] = rtt[tuple(rbt[j][pkk])]
    return (vbb,mbb,tuple(sbb),np.transpose(rbt))

# listVariablesListTransformRepasSort :: V.Vector Variable -> V.Vector TransformRepa -> V.Vector TransformRepa 

def listVariablesListTransformRepasSort(vv,ff):
    def der(tr):
        return set([transformRepasVarDerived(tr)])
    def vars(tr):
        return set(transformRepasVectorVar(tr))
    def und(tr):
        return vars(tr) - der(tr)
    def next(vv,ff):
        gg = []
        finished = False
        while not finished:
            found = False
            for (i,(tt,xx)) in enumerate(ff):
                if xx.issubset(vv):
                    found = True
                    break
            if found:
                vv = vv|der(tt)
                ff = ff[:i]+ff[i+1:]
                gg.append(tt)
            else:
                finished = True
        return gg
    vv1 = set(vv)
    ff1 = [(tt,und(tt)) for tt in ff if not der(tt).issubset(vv1)]
    return next(vv1,ff1)

# listVariablesListTransformRepasSort_1 :: V.Vector Variable -> V.Vector TransformRepa -> V.Vector TransformRepa 

def listVariablesListTransformRepasSort_1(vv,ff):
    def der(tr):
        return set([transformRepasVarDerived(tr)])
    def vars(tr):
        return set(transformRepasVectorVar(tr))
    def und(tr):
        return vars(tr) - der(tr)
    def next(vv,ff,gg):
        found = False
        for (i,(tt,xx)) in enumerate(ff):
            if xx.issubset(vv):
                found = True
                break
        if found:
            return next(vv|der(tt),ff[:i]+ff[i+1:],gg+[tt])
        return gg
    vv1 = set(vv)
    ff1 = [(tt,und(tt)) for tt in ff if not der(tt).issubset(vv1)]
    return next(vv1,ff1,[])

# historyRepasListTransformRepasApply :: HistoryRepa -> V.Vector TransformRepa -> HistoryRepa 

def historyRepasListTransformRepasApply(hr,fr):
    vars = historyRepasVectorVar
    ltrsort = listVariablesListTransformRepasSort
    ltrmul = historyRepasListTransformRepasApply_u
    return ltrmul(hr,ltrsort(vars(hr),fr))

# systemsFudsHistoryRepasMultiply_u :: System -> Fud -> HistoryRepa -> HistoryRepa 

def systemsFudsHistoryRepasMultiply_u(uu,ff,hr):
    ffqq = fudsSetTransform
    tttr = systemsTransformsTransformRepa_u
    ltrmul = historyRepasListTransformRepasApply
    gr = [tttr(uu,tt) for tt in ffqq(ff)]
    return ltrmul(hr,gr)

# systemsDecompFudsHistoryRepasMultiply :: System -> DecompFud -> HistoryRepa -> Tree ((State,Fud),HistoryRepa)

def systemsDecompFudsHistoryRepasMultiply(uu,df,aa):
    def unit(ss):
        return setStatesHistogramUnit(sset([ss]))
    fder = fudsDerived
    aahh = histogramsHistory
    hhhr = systemsHistoriesHistoryRepa
    empty = historyRepaEmpty
    size = historyRepasSize
    vars = historyRepasSetVariable
    def select(uu,ss,hh):
        return historyRepasHistoryRepasHistoryRepaSelection_u(hhhr(uu,aahh(unit(ss))),hh)
    def red(hr,vv):
        return setVarsHistoryRepasHistoryRepaReduced(vv,hr)
    def fmul(uu,ff,hh):
        return systemsFudsHistoryRepasMultiply_u(uu,ff,hh)
    def apply(zz,vv,aa):
        ll = []
        if size(aa) > 0:
            for ((ss,ff),yy) in zz.items():
                aa1 = select(uu,ss,aa)
                ww = fder(ff)
                bb = empty()
                if size(aa1) > 0:
                    bb = red(fmul(uu,ff,aa1),vv|ww)
                ll.append((((ss,ff),hlist(bb)),apply(yy,vv,bb)))
        return sdict(ll)
    return apply(df,vars(aa),aa)

# systemsDecompFudsHistoryRepasDecompFudReduced :: System -> DecompFud -> HistoryRepa -> DecompFud

def systemsDecompFudsHistoryRepasDecompFudReduced(uu,df,aa):
    fder = fudsDerived
    def red(ss,v):
        return setVarsStatesStateFiltered(sset([v]),ss)
    def fdep(ff,x):
        return fudsSetVarsDepends(ff,sset([x]))
    hrsize = historyRepasSize
    hrmult = systemsDecompFudsHistoryRepasMultiply
    def apply(w,zz):
        qq = sdict()
        for (((ss,ff),hr),yy) in zz.items():
            u = list(fder(ff))[0]
            xx = (red(ss,w),fdep(ff,u))
            a = hrsize(hr)
            if xx in qq:
                (b,_) = qq[xx]
                if a > b:
                    qq[xx] = (a,apply(u,yy))
            else:
                qq[xx] = (a,apply(u,yy))
        return sdict([(x,yy) for (x,(_,yy)) in qq.items()])
    (((ss,ff),_),yy) = list(hrmult(uu,df,aa).items())[0]
    w = list(fder(ff))[0]
    df1 = sdict([((ss,fdep(ff,w)),apply(w,yy))])
    return df1

# int listVarsArrayHistoriesAlignedTop_u(
#     int xmax, int omax, int n, int* svv, int m, int z1, int z2,
#     int* ppww, int* phh1, double* pxx1, int* phh2, double* pxx2,
#     int* tww1, int* tww2, double* ts1, double* ts2, int* ts3, int* s)

# parametersSetVarsHistoryRepasSetSetVarsAlignedTop_u :: Integer -> Integer -> Set.Set Variable -> HistoryRepa -> HistogramRepaRed -> HistoryRepa -> HistogramRepaRed -> (V.Vector ((Double,Double,Integer),Set.Set Variable),Integer) 

def parametersSetVarsHistoryRepasSetSetVarsAlignedTop_u(xmax,omax,ww,hh,hhx,hhrr,hhrrx):
    (vhh,mvv,svv,phh1) = hh
    (_,_,_,lxx1) = hhx
    (_,_,_,phh2) = hhrr
    (_,_,_,lxx2) = hhrrx
    (n,z) = phh1.shape
    (_,zrr) = phh2.shape
    vww = list(ww)
    m = len(vww)
    pww = [mvv[v] for v in vww]
    psvv = np.array(svv,dtype=np.dtype(np.int32))
    ppww = np.array(pww,dtype=np.dtype(np.int32))
    pxx1 = np.concatenate(lxx1)
    pxx2 = np.concatenate(lxx2)
    tww1 = np.zeros((omax),dtype=np.dtype(np.int32))
    tww2 = np.zeros((omax),dtype=np.dtype(np.int32))
    ts1 = np.zeros((omax),dtype=np.dtype(np.float64))
    ts2 = np.zeros((omax),dtype=np.dtype(np.float64))
    ts3 = np.zeros((omax),dtype=np.dtype(np.int32))
    s = np.zeros((1),dtype=np.dtype(np.int32))
    t = np.zeros((1),dtype=np.dtype(np.int32))
    listVarsArrayHistoriesAlignedTop_u_py(xmax,omax,n,psvv,m,z,zrr,ppww,phh1,pxx1,phh2,pxx2,tww1,tww2,ts1,ts2,ts3,s,t)
    t1 = t[0]
    qq = list(zip(list(zip(ts1[:t1],ts2[:t1],ts3[:t1])),[sset([vhh[p1],vhh[p2]]) for (p1,p2) in zip(tww1[:t1],tww2[:t1])]))
    return (qq,s[0])

if not AlignmentForeignPy_ok:
    parametersSetVarsHistoryRepasSetSetVarsAlignedTop_u = None


# int listVarsListTuplesArrayHistoriesAlignedTop_u(
#     int dense,
#     int xmax, int omax, int n, int* svv, int m, int d, int e,
#     int z1, int z2,
#     int* ppww, int* ppdd,
#     int* phh1, double* pxx1, int* phh2, double* pxx2,
#     int* tww1, int* tww2, double* ts1, double* ts2, int* ts3, int* s)

# parametersSetVarsSetSetVarsHistoryRepasSetSetVarsAlignedTop_u :: Integer -> Integer -> Set.Set Variable -> V.Vector (Set.Set Variable) -> HistoryRepa -> HistogramRepaRed -> HistoryRepa -> HistogramRepaRed -> (V.Vector ((Double,Double,Integer),Set.Set Variable),Integer) 

def parametersSetVarsSetSetVarsHistoryRepasSetSetVarsAlignedTop_u(xmax,omax,ww,vdd,hh,hhx,hhrr,hhrrx):
    (vhh,mvv,svv,phh1) = hh
    (_,_,_,lxx1) = hhx
    (_,_,_,phh2) = hhrr
    (_,_,_,lxx2) = hhrrx
    (n,z) = phh1.shape
    (_,zrr) = phh2.shape
    vww = list(ww)
    m = len(vww)
    d = len(vdd)
    e = len(vdd[0])
    pww = [mvv[v] for v in vww]
    pdd = [mvv[v] for qq in vdd for v in qq]
    psvv = np.array(svv,dtype=np.dtype(np.int32))
    ppww = np.array(pww,dtype=np.dtype(np.int32))
    ppdd = np.array(pdd,dtype=np.dtype(np.int32))
    pxx1 = np.concatenate(lxx1)
    pxx2 = np.concatenate(lxx2)
    tww1 = np.zeros((omax),dtype=np.dtype(np.int32))
    tww2 = np.zeros((omax),dtype=np.dtype(np.int32))
    ts1 = np.zeros((omax),dtype=np.dtype(np.float64))
    ts2 = np.zeros((omax),dtype=np.dtype(np.float64))
    ts3 = np.zeros((omax),dtype=np.dtype(np.int32))
    s = np.zeros((1),dtype=np.dtype(np.int32))
    t = np.zeros((1),dtype=np.dtype(np.int32))
#   listVars          ArrayHistoriesAlignedTop_u_py(  xmax,omax,n,psvv,m,    z,zrr,ppww,     phh1,pxx1,phh2,pxx2,tww1,tww2,ts1,ts2,ts3,s,t)
    listVarsListTuplesArrayHistoriesAlignedTop_u_py(0,xmax,omax,n,psvv,m,d,e,z,zrr,ppww,ppdd,phh1,pxx1,phh2,pxx2,tww1,tww2,ts1,ts2,ts3,s,t) 
    t1 = t[0]
    qq = list(zip(list(zip(ts1[:t1],ts2[:t1],ts3[:t1])),[vdd[p2]|sset([vww[p1]]) for (p1,p2) in zip(tww1[:t1],tww2[:t1])]))
    return (qq,s[0])

if not AlignmentForeignPy_ok:
    parametersSetVarsSetSetVarsHistoryRepasSetSetVarsAlignedTop_u = None


# int listVarsListTuplesArrayHistoriesAlignedExcludeHiddenTop_u(
#     int dense,
#     int xmax, int omax, int n, int* svv, int m, int d, int e,
#     int z1, int z2,
#     int ccl, int* ppccd, int* ppccu,
#     int* ppww, int* ppdd,
#     int* phh1, double* pxx1, int* phh2, double* pxx2,
#     int* tww1, int* tww2, double* ts1, double* ts2, int* ts3, int* s)

# parametersSetVarsSetSetVarsHistoryRepasSetSetVarsAlignedExcludeHiddenDenseTop_u :: Integer -> Integer -> Set.Set (Variable,Variable) -> Set.Set Variable -> V.Vector (Set.Set Variable) -> HistoryRepa -> HistogramRepaRed -> HistoryRepa -> HistogramRepaRed -> (V.Vector ((Double,Double,Integer),Set.Set Variable),Integer) 

def parametersSetVarsSetSetVarsHistoryRepasSetSetVarsAlignedExcludeHiddenDenseTop_u(wmax,omax,cc,ww,vdd,hh,hhx,hhrr,hhrrx):
    (vhh,mvv,svv,phh1) = hh
    (_,_,_,lxx1) = hhx
    (_,_,_,phh2) = hhrr
    (_,_,_,lxx2) = hhrrx
    (n,z) = phh1.shape
    (_,zrr) = phh2.shape
    vcc = list(cc)
    ccl = len(vcc)
    pccd = [mvv[d] for (d,u) in vcc]
    pccu = [mvv[u] for (d,u) in vcc]
    vww = list(ww)
    m = len(vww)
    d = len(vdd)
    e = len(vdd[0])
    pww = [mvv[v] for v in vww]
    pdd = [mvv[v] for qq in vdd for v in qq]
    psvv = np.array(svv,dtype=np.dtype(np.int32))
    ppccd = np.array(pccd,dtype=np.dtype(np.int32))
    ppccu = np.array(pccu,dtype=np.dtype(np.int32))
    ppww = np.array(pww,dtype=np.dtype(np.int32))
    ppdd = np.array(pdd,dtype=np.dtype(np.int32))
    pxx1 = np.concatenate(lxx1)
    pxx2 = np.concatenate(lxx2)
    tww1 = np.zeros((omax),dtype=np.dtype(np.int32))
    tww2 = np.zeros((omax),dtype=np.dtype(np.int32))
    ts1 = np.zeros((omax),dtype=np.dtype(np.float64))
    ts2 = np.zeros((omax),dtype=np.dtype(np.float64))
    ts3 = np.zeros((omax),dtype=np.dtype(np.int32))
    s = np.zeros((1),dtype=np.dtype(np.int32))
    t = np.zeros((1),dtype=np.dtype(np.int32))
#   listVars          ArrayHistoriesAligned             Top_u_py(  xmax,omax,n,psvv,m,    z,zrr,                ppww,     phh1,pxx1,phh2,pxx2,tww1,tww2,ts1,ts2,ts3,s,t)
#   listVarsListTuplesArrayHistoriesAligned             Top_u_py(0,xmax,omax,n,psvv,m,d,e,z,zrr,                ppww,ppdd,phh1,pxx1,phh2,pxx2,tww1,tww2,ts1,ts2,ts3,s,t) 
    listVarsListTuplesArrayHistoriesAlignedExcludeHiddenTop_u_py(1,wmax,omax,n,psvv,m,d,e,z,zrr,ccl,ppccd,ppccu,ppww,ppdd,phh1,pxx1,phh2,pxx2,tww1,tww2,ts1,ts2,ts3,s,t) 
    t1 = t[0]
    qq = list(zip(list(zip(ts1[:t1],ts2[:t1],ts3[:t1])),[vdd[p2]|sset([vww[p1]]) for (p1,p2) in zip(tww1[:t1],tww2[:t1])]))
    return (qq,s[0])

if not AlignmentForeignPy_ok:
    parametersSetVarsSetSetVarsHistoryRepasSetSetVarsAlignedTop_u = None


# void listListVarsArrayHistoryPairsPartitionIndependent_u(
#     double z, int v, int n, int* svv, int m, int r,
#     int* lyy, int* syy, int* pppp, double* aa1, double* aa2,
#     double* bb1, double* bb2)

# setSetVarsHistogramRepaPairPartitionIndependentPair_u :: 
#   Set.Set (Set.Set Variable) -> HistogramRepaVec -> HistogramRepaVec 

def setSetVarsHistogramRepaPairPartitionIndependentPair_u(pp,rrv):
    def prod(ll):
        p = 1
        for x in ll:
            p *= x
        return p
    (vvv,mvv,z,svv,vaa) = rrv
    [aa,aarr] = vaa
    v = prod(svv)
    n = len(svv)
    vpp = [list(cc) for cc in pp]
    m = len(vpp)
    vyy = [VarIndex(i) for i in range(m)]
    myy = sdict([(v,i) for (i,v) in enumerate(vyy)])
    ppp = [[mvv[v] for v in vcc] for vcc in vpp]
    ppp1 = [p for pcc in ppp for p in pcc]
    syy = tuple([prod([svv[p] for p in pcc]) for pcc in ppp])
    lyy = [len(pcc) for pcc in ppp]
    r = sum(syy)
    psvv = np.array(svv,dtype=np.dtype(np.int32))
    plyy = np.array(lyy,dtype=np.dtype(np.int32))
    psyy = np.array(syy,dtype=np.dtype(np.int32))
    pppp = np.array(ppp1,dtype=np.dtype(np.int32))
    bb = np.zeros(syy,dtype=np.dtype(np.float64))
    bbrr = np.zeros(syy,dtype=np.dtype(np.float64))
    listListVarsArrayHistoryPairsPartitionIndependent_u_py(z,v,n,psvv,m,r,plyy,psyy,pppp,aa,aarr,bb,bbrr)
    vbb = [bb,bbrr]
    return (vyy,myy,z,syy,vbb)

if not AlignmentForeignPy_ok:
    setSetVarsHistogramRepaPairPartitionIndependentPair_u = None


# int listListVarsArrayHistoryPairsSetTuplePartitionTop_u(
#     int pmax, double z, int v, int n, int* svv, int q, double y1,
#     int* qm, int* ql, int* qs, int* qp, double* aa1, double* aa2,
#     int* tt)

# parametersHistogramRepaVecsSetTuplePartitionTop_u :: 
#   Integer -> Integer -> Integer -> HistogramRepaVec -> Double -> (Set.Set (Set.Set (Set.Set Variable)),Integer)

def parametersHistogramRepaVecsSetTuplePartitionTop_u(mmax,umax,pmax,rrv,y1):
    def prod(ll):
        p = 1
        for x in ll:
            p *= x
        return p
    def nlluv(n,ll):
        d = n - len(ll)
        if d > 0:
            return ll + [0]*d
        else:
            return ll[:n]
    (vvv,_,z,svv,vaa) = rrv
    [aa, _, aarr, _] = vaa
    v = prod(svv)
    n = len(svv)
    mm = [[0]]
    for i in range(2,n+1):
        mm1 = []
        for xx in mm:
            for j in range(i):
                if j < mmax and j <= max(xx) + 1:
                    mm1.append([j]+xx)
        mm = mm1
    qq = []
    for ll in mm[1:]:
        pp = sdict()
        for (p,c) in enumerate(ll):
            if c not in pp:
                pp[c] = sset([p])
            pp[c].add(p)
        pp1 = [list(cc) for (_,cc) in pp.items()]
        rr = [prod([svv[p] for p in pcc]) for pcc in pp1]
        if all([u <= umax for u in rr]):
            qq.append((pp1,rr))
    q = len(qq)
    qm = [len(pp) for (pp,_) in qq]
    ql = [i for (pp,_) in qq for i in nlluv(n,[len(cc) for cc in pp])]
    qs = [i for (_,rr) in qq for i in nlluv(n,rr)]
    qp = [i for (pp,_) in qq for i in nlluv(n,[j for cc in pp for j in cc])]
    psvv = np.array(svv,dtype=np.dtype(np.int32))
    pqm = np.array(qm,dtype=np.dtype(np.int32))
    pql = np.array(ql,dtype=np.dtype(np.int32))
    pqs = np.array(qs,dtype=np.dtype(np.int32))
    pqp = np.array(qp,dtype=np.dtype(np.int32))
    ptt = np.zeros((pmax),dtype=np.dtype(np.int32))
    t = np.zeros((1),dtype=np.dtype(np.int32))
    listListVarsArrayHistoryPairsSetTuplePartitionTop_u_py(pmax,z,v,n,psvv,q,y1,pqm,pql,pqs,pqp,aa,aarr,ptt,t)
    t1 = t[0] 
    tt1 = sset([sset([sset([vvv[j] for j in cc]) for cc in qq[i][0]]) for i in ptt[:t1]])
    return (tt1,q)

if not AlignmentForeignPy_ok:
    parametersHistogramRepaVecsSetTuplePartitionTop_u = None

# int listListVarsArrayHistoryPairsSetTuplePartitionTop_u(
#     int pmax, double z, int v, int n, int* svv, int q, double y1,
#     int* qm, int* ql, int* qs, int* qp, double* aa1, double* aa2,
#     int* tt)

# parametersHistogramRepaVecsSetTuplePartitionTopByM_u :: 
#   Integer -> Integer -> Integer -> HistogramRepaVec -> Double -> (Set.Set (Set.Set (Set.Set Variable)),Integer)

def parametersHistogramRepaVecsSetTuplePartitionTopByM_u(mmax,umax,pmax,rrv,y1):
    def prod(ll):
        p = 1
        for x in ll:
            p *= x
        return p
    def nlluv(n,ll):
        d = n - len(ll)
        if d > 0:
            return ll + [0]*d
        else:
            return ll[:n]
    (vvv,_,z,svv,vaa) = rrv
    [aa, _, aarr, _] = vaa
    v = prod(svv)
    n = len(svv)
    def parter(qq):
        q = len(qq)
        qm = [len(pp) for (pp,_) in qq]
        ql = [i for (pp,_) in qq for i in nlluv(n,[len(cc) for cc in pp])]
        qs = [i for (_,rr) in qq for i in nlluv(n,rr)]
        qp = [i for (pp,_) in qq for i in nlluv(n,[j for cc in pp for j in cc])]
        psvv = np.array(svv,dtype=np.dtype(np.int32))
        pqm = np.array(qm,dtype=np.dtype(np.int32))
        pql = np.array(ql,dtype=np.dtype(np.int32))
        pqs = np.array(qs,dtype=np.dtype(np.int32))
        pqp = np.array(qp,dtype=np.dtype(np.int32))
        ptt = np.zeros((pmax),dtype=np.dtype(np.int32))
        t = np.zeros((1),dtype=np.dtype(np.int32))
        listListVarsArrayHistoryPairsSetTuplePartitionTop_u_py(pmax,z,v,n,psvv,q,y1,pqm,pql,pqs,pqp,aa,aarr,ptt,t)
        t1 = t[0] 
        return sset([sset([sset([vvv[j] for j in cc]) for cc in qq[i][0]]) for i in ptt[:t1]])
    mm = [[0]]
    for i in range(2,n+1):
        mm1 = []
        for xx in mm:
            for j in range(i):
                if j < mmax and j <= max(xx) + 1:
                    mm1.append([j]+xx)
        mm = mm1
    qq = []
    for ll in mm[1:]:
        pp = sdict()
        for (p,c) in enumerate(ll):
            if c not in pp:
                pp[c] = sset([p])
            pp[c].add(p)
        pp1 = [list(cc) for (_,cc) in pp.items()]
        rr = [prod([svv[p] for p in pcc]) for pcc in pp1]
        if all([u <= umax for u in rr]):
            qq.append((pp1,rr))
    tt1 = sset()
    for m in range(2,mmax+1):
        qq1 = [(pp,rr) for (pp,rr) in qq if len(pp) == m]
        tt1 |= parter(qq1)
    return (tt1,len(qq))

if not AlignmentForeignPy_ok:
    parametersHistogramRepaVecsSetTuplePartitionTopByM_u = None


# int arrayHistoryPairsRollMax_u(
#     int v, int n, int* svvy, int d, int nd,
#     double* aay, double* aaxy, double* bby, double* bbxy,
#     int* ppm)

# histogramRepaVecsRollMax :: HistogramRepaVec -> (V.Vector (UV.Vector Int),Integer)

def histogramRepaVecsRollMax(rrv):
    def prod(ll):
        p = 1
        for x in ll:
            p *= x
        return p
    (_,_,_,svv,vaa) = rrv
    [aa,aax,bb,bbx] = vaa
    v = prod(svv)
    n = len(svv)
    d = max(svv)
    nd = n * d
    psvv = np.array(svv,dtype=np.dtype(np.int32))
    ppm = np.zeros((nd,),dtype=np.dtype(np.int32))
    t = np.zeros((1),dtype=np.dtype(np.int32))
    arrayHistoryPairsRollMax_u_py(v,n,psvv,d,nd,aa,aax,bb,bbx,ppm,t)
    t1 = t[0] 
    tt = [list(ppm[d*i:d*i+e]) for (i,e) in enumerate(svv)]
    return (tt,t1)

if not AlignmentForeignPy_ok:
    histogramRepaVecsRollMax = None
