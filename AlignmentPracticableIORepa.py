from AlignmentPracticableRepa import *
import logging
from timeit import default_timer as timer
from sys import stdout

# logging.basicConfig(format='%(asctime)s : %(name)s : %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=logging.INFO)
# logging.basicConfig(format='%(message)s', level=logging.INFO)
logging.basicConfig(format='%(message)s')

layerer_log = logging.getLogger('layerer')
layerer_log.setLevel(logging.INFO)
tupler_log = logging.getLogger('tupler')
tupler_log.setLevel(logging.INFO)
parter_log = logging.getLogger('parter')
parter_log.setLevel(logging.INFO)
roller_log = logging.getLogger('roller')
roller_log.setLevel(logging.INFO)
applier_log = logging.getLogger('applier')
applier_log.setLevel(logging.INFO)
dervarser_log = logging.getLogger('dervarser')
dervarser_log.setLevel(logging.INFO)

decomper_log = logging.getLogger('decomper')
decomper_log.setLevel(logging.INFO)


# parametersSystemsLayererMaxRollByMExcludedSelfHighestIORepa_u :: 
#   Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> 
#   System -> Set.Set Variable -> HistoryRepa -> HistogramRepaRed -> HistoryRepa -> HistogramRepaRed -> Integer ->
#   IO (System, Fud, [(Set.Set Variable, Double)])

def parametersSystemsLayererMaxRollByMExcludedSelfHighestIORepa_u(wmax,lmax,xmax,omax,bmax,mmax,umax,pmax,uu,vv,xx,xxp,xxrr,xxrrp,f):
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
        return parametersSystemsBuilderTupleNoSumlayerMultiEffectiveRepa_ui(xmax,omax,bmax,mmax,uu,vv,ff,hh,hhp,hhrr,hhrrp)
    def parter(uu,kk,bb,y1):
        return parametersSystemsPartitionerMaxRollByMRepa_ui(mmax,umax,pmax,uu,kk,bb,y1)
    def roller(qq):
        return parametersRollerMaximumRollExcludedSelfRepa_i(qq)
    def buildffdervar(uu,vv,ff,xx,xxp,xxrr,xxrrp):
        (x1,s1) = parametersSystemsBuilderDerivedVarsHighestNoSumlayerRepa_ui(wmax,omax,uu,vv,ff,xx,xxp,xxrr,xxrrp)
        return ([(kk,a) for ((kk,_,_),a) in x1],s1)
    def layer(vv,uu,ff,mm,xx,xxp,xxrr,xxrrp,f,l):
        if l > lmax:
            return (uu,ff,mm)
        layerer_log.info(">>> layer\tfud: %d\tlayer: %d" % (f,l))
        t1 = timer()
        tupler_log.info(">>> tupler")
        tupler_log.info("substrate cardinality: %d" % len(vv))
        tupler_log.info("fud cardinality: %d" % len(ffqq(ff)))
        stdout.flush()
        (x2,s2) = buildfftup(uu,vv,ff,xx,xxp,xxrr,xxrrp)
        if len(x2) > 0:
            tupler_log.info("tuple cardinality: %d" % len(x2))
            tupler_log.info("max tuple algn: %.2f" % max([b for (a,b) in x2]))
        else:
            tupler_log.info("no tuples")
        t2 = timer()
        tupler_log.info("tupler\tsearched: %d\trate: %.2f" % (s2,s2/(t2-t1)))
        tupler_log.info("<<< tupler %.3fs" % (t2-t1))
        parter_log.info(">>> parter")
        stdout.flush()
        y3 = [parter(uu,kk,bb,y1) for ((kk,bb),y1) in x2]
        x3 = [x for (ll,_) in y3 for x in ll]
        s3 = sum([s for (_,s) in y3])
        if len(x3) > 0:
            parter_log.info("partitions cardinality: %d" % len(x3))
        else:
            parter_log.info("no tuple partitions")
        t3 = timer()
        parter_log.info("parter\tsearched: %d\trate: %.2f" % (s3,s3/(t3-t2)))
        parter_log.info("<<< parter %.3fs" % (t3-t2))
        roller_log.info(">>> roller")
        stdout.flush()
        y4 = [roller(qq) for qq in x3]
        x4 = [x for (ll,_) in y4 for x in ll]
        s4 = sum([s for (_,s) in y4])
        if len(x4) > 0:
            roller_log.info("roll cardinality: %d" % len(x4))
        else:
            roller_log.info("no rolls")
        t4 = timer()
        roller_log.info("roller\tsearched: %d\trate: %.2f" % (s4,s4/(t4-t3)))
        roller_log.info("<<< roller %.3fs" % (t4-t3))
        applier_log.info(">>> application")
        stdout.flush()
        ll0 = []
        for (yy,pp) in x4:
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
            applier_log.info("fud cardinality: %d" % len(ffqq(gg)))
            t5 = timer()
            applier_log.info("<<< application %.3fs" % (t5-t4))
            dervarser_log.info( ">>> dervarser")
            stdout.flush()
            (mm1,s5) = buildffdervar(uu1,vv,gg,xx1,xxp1,xxrr1,xxrrp1)
            if len(mm1) > 0:
                dervarser_log.info("der vars algn density: %.2f" % maxr(mm1))
            else:
                dervarser_log.info("no der vars sets")
            t6 = timer()
            dervarser_log.info("dervarser\tsearched: %d\trate: %.2f" % (s5,s5/(t6-t5)))
            dervarser_log.info("<<< dervarser %.3fs" % (t6-t5))
            layerer_log.info( "<<< layer %.3fs" % (t6-t1))
            stdout.flush()
            if l <= lmax and (len(mm) == 0 or maxr(mm1) > maxr(mm) + repaRounding):
                (ffr,ll0,ll,ll1) = (None,None,None,None)
                (x2,x3,x4) = (None,None,None)
                return layer(vv,uu1,gg,mm1,xx1,xxp1,xxrr1,xxrrp1,f,l+1)
        else:
            t5 = timer()
            applier_log.info("<<< application %.3fs" % (t5-t4))
            layerer_log.info( "<<< layer %.3fs" % (t5-t1))
            stdout.flush()
        return (uu,ff,mm) 
    layerer_log.info(">>> layerer")
    t1 = timer()
    x1 = layer(vv,uu,fudEmpty(),[],xx,xxp,xxrr,xxrrp,f,1)
    t2 = timer()
    layerer_log.info("<<< layerer %.3fs" % (t2-t1))
    stdout.flush()
    return x1

# parametersSystemsHistoryRepasDecomperMaxRollByMExcludedSelfHighestFmaxIORepa :: 
#   Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> 
#   Integer -> Integer ->
#   System -> Set.Set Variable -> HistoryRepa -> 
#   IO (Maybe (System, DecompFud))

def parametersSystemsHistoryRepasDecomperMaxRollByMExcludedSelfHighestFmaxIORepa(wmax,lmax,xmax,omax,bmax,mmax,umax,pmax,fmax,mult,seed,uu,vv,aa):
    repaRounding = 1e-6
    dom = relationsDomain
    def maxd(mm):
        if len(mm) > 0:
            return list(sset([(b,a) for (a,b) in mm]))[-1]
        return (0,sset())
    def tsgl(r):
        return sdict([(r,sdict())])
    uvars = systemsSetVar
    acard = histogramsCardinality
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
    dfzz = decompFudsTreePairStateFud
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
        decomper_log.info(">>> repa shuffle")
        stdout.flush()
        t1 = timer()
        z = size(xx)
        xxrr = hrconcat([hrshuffle(xx,seed+i*z) for i in range(1,mult+1)])
        t2 = timer()
        decomper_log.info("<<< repa shuffle %.3fs" % (t2-t1))
        decomper_log.info(">>> repa perimeters")
        stdout.flush()
        t1 = timer()
        xxp = hrhx(xx)
        xxrrp = hrhx(xxrr)
        t2 = timer()
        decomper_log.info("<<< repa perimeters %.3fs" % (t2-t1))
        return parametersSystemsLayererMaxRollByMExcludedSelfHighestIORepa_u(wmax,lmax,xmax,omax,bmax,mmax,umax,pmax,uu,vv,xx,xxp,xxrr,xxrrp,f)
    def decomp(uu,zz,qq,f):
        if len(zz) == 0:
            (uur,ffr,nnr) = layerer(uu,aa,f)
            if len(ffr) == 0 or len(nnr) == 0:
                return (uu, decompFudEmpty())
            (ar,kkr) = maxd(nnr)
            if ar <= repaRounding:
                return (uu, decompFudEmpty())
            decomper_log.info(">>> slicing")
            stdout.flush()
            t3 = timer()
            ffr1 = depends(ffr,kkr)
            decomper_log.info("dependent fud cardinality : %d" % len(ffqq(ffr1)))
            aar = apply(uur,ffr1,aa)
            aa1 = trim(reduce(uur,fder(ffr1),aar))
            decomper_log.info("derived cardinality : %d" % acard(aa1))
            zzr = tsgl((stateEmpty(),ffr1))
            qq[(stateEmpty(),ffr1)] = (aar,aa1)
            (ffr,nnr,kkr) = (None,None,None)
            t4 = timer()
            decomper_log.info("<<< slicing %.3fs" % (t4-t3))
            stdout.flush()
            return decomp(uur,zzr,qq,f+1)
        if fmax > 0 and f > fmax:
            return (uu,zzdf(zztrim(zz)))
        decomper_log.info(">>> slice selection")
        stdout.flush()
        t1 = timer()
        mm = []
        for (nn,yy) in treesPlaces(zz):
            (rr,ff) = nn[-1]
            if len(ff) > 0:
                (bb,bb1) = qq[(rr,ff)]
                tt = dom(treesRoots(yy))
                for (ss,a) in aall(red(bb1,fder(ff))):
                    if a > 0 and ss not in tt:
                        mm.append((a,(nn,ss,bb)))
        decomper_log.info("slices: %d" % len(mm))
        if len(mm) == 0:
            t2 = timer()
            decomper_log.info("<<< slice selection %.3fs" % (t2-t1))
            stdout.flush()
            return (uu,zzdf(zztrim(zz)))
        mm.sort(key = lambda x: x[0])
        (a,(nn,ss,bb)) = mm[-1]
        cc = hrhrred(select(uu,ss,bb),vars(aa))
        decomper_log.info("decomp path length: %d" % len(nn))
        decomper_log.info("slice size: %d" % a)
        t2 = timer()
        decomper_log.info("<<< slice selection %.3fs" % (t2-t1))
        stdout.flush()
        (uuc,ffc,nnc) = layerer(uu,cc,f)
        decomper_log.info(">>> slicing")
        stdout.flush()
        t3 = timer()
        (ac,kkc) = maxd(nnc)
        ffc1 = fudEmpty()
        if ac > repaRounding:
             ffc1 = depends(ffc,kkc)
        decomper_log.info("dependent fud cardinality : %d" % len(ffqq(ffc1)))
        ccc = apply(uuc,ffc1,cc)
        cc1 = trim(reduce(uuc,fder(ffc1),ccc))
        decomper_log.info("derived cardinality : %d" % acard(cc1))
        qq[(ss,ffc1)] = (ccc,cc1)
        zzc = pathsTree(treesPaths(zz) + [nn+[(ss,ffc1)]])
        (mm,cc,ffc,nnc,kkc) = (None,None,None,None,None)
        t4 = timer()
        decomper_log.info("<<< slicing %.3fs" % (t4-t3))
        stdout.flush()
        return decomp(uuc,zzc,qq,f+1)
    if wmax < 0 or lmax < 0 or xmax < 0 or omax < 0 or bmax < 0 or mmax < 1 or umax < 0 or pmax < 0:
        return None
    if size(aa) == 0 or mult < 1:
        return None
    if not (vars(aa).issubset(uvars(uu)) and vv.issubset(vars(aa))):
        return None
    decomper_log.info(">>> decomper")
    t1 = timer()
    x1 = decomp(uu,emptyTree(),sdict(),1)
    decomper_log.info("nodes: %d" % len(treesNodes(dfzz(x1[1]))))
    t2 = timer()
    decomper_log.info("<<< decomper repa %.3fs" % (t2 - t1))
    stdout.flush()
    return x1

# parametersSystemsLayererLevelMaxRollByMExcludedSelfHighestIORepa_u :: 
#   Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> 
#   System -> Set.Set Variable -> Fud -> 
#   HistoryRepa -> HistogramRepaRed -> HistoryRepa -> HistogramRepaRed -> Integer -> Integer ->
#   IO (System, Fud, [(Set.Set Variable, Double)])

def parametersSystemsLayererLevelMaxRollByMExcludedSelfHighestIORepa_u(wmax,lmax,xmax,omax,bmax,mmax,umax,pmax,uu,vvg,ffg,xx,xxp,xxrr,xxrrp,f,g):
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
    fund = fudsUnderlying
    fder = fudsDerived
    depends = fudsSetVarsDepends
    def funion(ff,gg):
        return qqff(ffqq(ff) | ffqq(gg))
    def buildfftup(uu,vvg,ffg,ff,hh,hhp,hhrr,hhrrp):
        return parametersSystemsBuilderTupleLevelNoSumlayerMultiEffectiveRepa_ui(xmax,omax,bmax,mmax,uu,vvg,ffg,ff,hh,hhp,hhrr,hhrrp)
    def parter(uu,kk,bb,y1):
        return parametersSystemsPartitionerMaxRollByMRepa_ui(mmax,umax,pmax,uu,kk,bb,y1)
    def roller(qq):
        return parametersRollerMaximumRollExcludedSelfRepa_i(qq)
    def buildffdervar(uu,vv,ffg,ff,xx,xxp,xxrr,xxrrp):
        (x1,s1) = parametersSystemsBuilderDerivedVarsLevelHighestNoSumlayerRepa_ui(wmax,omax,uu,vv,ffg,ff,xx,xxp,xxrr,xxrrp)
        return ([(kk,a) for ((kk,_,_),a) in x1],s1)
    def layer(uu,ff,mm,xx,xxp,xxrr,xxrrp,l):
        if l > lmax:
            return (uu,ff,mm)
        layerer_log.info(">>> layer\tfud: %d\tlevel node: %d\tlayer: %d" % (f,g,l))
        t1 = timer()
        tupler_log.info(">>> tupler")
        tupler_log.info("level substrate cardinality: %d" % len(vvg))
        tupler_log.info("level fud derived cardinality: %d" % len(fder(ffg)))
        tupler_log.info("fud cardinality: %d" % len(ffqq(ff)))
        tupler_log.info("level excluded fud cardinality: %d" % len(ffqq(ff)-ffqq(ffg)))
        stdout.flush()
        (x2,s2) = buildfftup(uu,vvg,ffg,ff,xx,xxp,xxrr,xxrrp)
        if len(x2) > 0:
            tupler_log.info("tuple cardinality: %d" % len(x2))
            tupler_log.info("max tuple algn: %.2f" % max([b for (a,b) in x2]))
        else:
            tupler_log.info("no tuples")
        t2 = timer()
        tupler_log.info("tupler\tsearched: %d\trate: %.2f" % (s2,s2/(t2-t1)))
        tupler_log.info("<<< tupler %.3fs" % (t2-t1))
        parter_log.info(">>> parter")
        stdout.flush()
        y3 = [parter(uu,kk,bb,y1) for ((kk,bb),y1) in x2]
        x3 = [x for (ll,_) in y3 for x in ll]
        s3 = sum([s for (_,s) in y3])
        if len(x3) > 0:
            parter_log.info("partitions cardinality: %d" % len(x3))
        else:
            parter_log.info("no tuple partitions")
        t3 = timer()
        parter_log.info("parter\tsearched: %d\trate: %.2f" % (s3,s3/(t3-t2)))
        parter_log.info("<<< parter %.3fs" % (t3-t2))
        roller_log.info(">>> roller")
        stdout.flush()
        y4 = [roller(qq) for qq in x3]
        x4 = [x for (ll,_) in y4 for x in ll]
        s4 = sum([s for (_,s) in y4])
        if len(x4) > 0:
            roller_log.info("roll cardinality: %d" % len(x4))
        else:
            roller_log.info("no rolls")
        t4 = timer()
        roller_log.info("roller\tsearched: %d\trate: %.2f" % (s4,s4/(t4-t3)))
        roller_log.info("<<< roller %.3fs" % (t4-t3))
        applier_log.info(">>> application")
        stdout.flush()
        ll0 = []
        for (yy,pp) in x4:
            for (jj,p) in zip(yy,pp):
                if max(p) + 1 < len(p):
                    ii = list(zip(cart(uu,jj),p))
                    ll0.append(ii)
        ll = []
        for (b,ii) in enumerate(ll0):
            w = VarPair((VarPair((VarPair((VarInt(f),VarInt(g))),VarInt(l))),VarInt(b+1)))
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
            gg = funion(funion(ff,hh),depends(ffg,fund(hh)))
            applier_log.info("fud cardinality: %d" % len(ffqq(gg)))
            t5 = timer()
            applier_log.info("<<< application %.3fs" % (t5-t4))
            dervarser_log.info( ">>> dervarser")
            stdout.flush()
            (mm1,s5) = buildffdervar(uu1,vvg,ffg,gg,xx1,xxp1,xxrr1,xxrrp1)
            if len(mm1) > 0:
                dervarser_log.info("der vars algn density: %.2f" % maxr(mm1))
            else:
                dervarser_log.info("no der vars sets")
            t6 = timer()
            dervarser_log.info("dervarser\tsearched: %d\trate: %.2f" % (s5,s5/(t6-t5)))
            dervarser_log.info("<<< dervarser %.3fs" % (t6-t5))
            layerer_log.info( "<<< layer %.3fs" % (t6-t1))
            stdout.flush()
            if l <= lmax and (len(mm) == 0 or maxr(mm1) > maxr(mm) + repaRounding):
                (ffr,ll0,ll,ll1) = (None,None,None,None)
                (x2,x3,x4) = (None,None,None)
                return layer(uu1,gg,mm1,xx1,xxp1,xxrr1,xxrrp1,l+1)
        else:
            t5 = timer()
            applier_log.info("<<< application %.3fs" % (t5-t4))
            layerer_log.info( "<<< layer %.3fs" % (t5-t1))
            stdout.flush()
        return (uu,ff,mm) 
    layerer_log.info(">>> layerer")
    t1 = timer()
    x1 = layer(uu,fudEmpty(),[],xx,xxp,xxrr,xxrrp,1)
    t2 = timer()
    layerer_log.info("<<< layerer %.3fs" % (t2-t1))
    stdout.flush()
    return x1

# TBD

# parametersSystemsHistoryRepasDecomperLevelMaxRollByMExcludedSelfHighestFmaxIORepa :: 
#   Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> Integer -> 
#   Integer -> Integer ->
#   System -> Tree (Integer, Set.Set Variable, Fud) -> HistoryRepa -> 
#   IO (Maybe (System, DecompFud))

def parametersSystemsHistoryRepasDecomperLevelMaxRollByMExcludedSelfHighestFmaxIORepa(wmax,lmax,xmax,omax,bmax,mmax,umax,pmax,fmax,mult,seed,uu,zzg,aa):
    repaRounding = 1e-6
    dom = relationsDomain
    def maxd(mm):
        if len(mm) > 0:
            return list(sset([(b,a) for (a,b) in mm]))[-1]
        return (0,sset())
    def tsgl(r):
        return sdict([(r,sdict())])
    uvars = systemsSetVar
    acard = histogramsCardinality
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
    dfzz = decompFudsTreePairStateFud
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
        decomper_log.info(">>> repa shuffle")
        stdout.flush()
        t1 = timer()
        z = size(xx)
        xxrr = hrconcat([hrshuffle(xx,seed+i*z) for i in range(1,mult+1)])
        t2 = timer()
        decomper_log.info("<<< repa shuffle %.3fs" % (t2-t1))
        decomper_log.info(">>> repa perimeters")
        stdout.flush()
        t1 = timer()
        xxp = hrhx(xx)
        xxrrp = hrhx(xxrr)
        t2 = timer()
        decomper_log.info("<<< repa perimeters %.3fs" % (t2-t1))
        return parametersSystemsLayererMaxRollByMExcludedSelfHighestIORepa_u(wmax,lmax,xmax,omax,bmax,mmax,umax,pmax,uu,vv,xx,xxp,xxrr,xxrrp,f)
    def decomp(uu,zz,qq,f):
        if len(zz) == 0:
            (uur,ffr,nnr) = layerer(uu,aa,f)
            if len(ffr) == 0 or len(nnr) == 0:
                return (uu, decompFudEmpty())
            (ar,kkr) = maxd(nnr)
            if ar <= repaRounding:
                return (uu, decompFudEmpty())
            decomper_log.info(">>> slicing")
            stdout.flush()
            t3 = timer()
            ffr1 = depends(ffr,kkr)
            decomper_log.info("dependent fud cardinality : %d" % len(ffqq(ffr1)))
            aar = apply(uur,ffr1,aa)
            aa1 = trim(reduce(uur,fder(ffr1),aar))
            decomper_log.info("derived cardinality : %d" % acard(aa1))
            zzr = tsgl((stateEmpty(),ffr1))
            qq[(stateEmpty(),ffr1)] = (aar,aa1)
            (ffr,nnr,kkr) = (None,None,None)
            t4 = timer()
            decomper_log.info("<<< slicing %.3fs" % (t4-t3))
            stdout.flush()
            return decomp(uur,zzr,qq,f+1)
        if fmax > 0 and f > fmax:
            return (uu,zzdf(zztrim(zz)))
        decomper_log.info(">>> slice selection")
        stdout.flush()
        t1 = timer()
        mm = []
        for (nn,yy) in treesPlaces(zz):
            (rr,ff) = nn[-1]
            if len(ff) > 0:
                (bb,bb1) = qq[(rr,ff)]
                tt = dom(treesRoots(yy))
                for (ss,a) in aall(red(bb1,fder(ff))):
                    if a > 0 and ss not in tt:
                        mm.append((a,(nn,ss,bb)))
        decomper_log.info("slices: %d" % len(mm))
        if len(mm) == 0:
            t2 = timer()
            decomper_log.info("<<< slice selection %.3fs" % (t2-t1))
            stdout.flush()
            return (uu,zzdf(zztrim(zz)))
        mm.sort(key = lambda x: x[0])
        (a,(nn,ss,bb)) = mm[-1]
        cc = hrhrred(select(uu,ss,bb),vars(aa))
        decomper_log.info("decomp path length: %d" % len(nn))
        decomper_log.info("slice size: %d" % a)
        t2 = timer()
        decomper_log.info("<<< slice selection %.3fs" % (t2-t1))
        stdout.flush()
        (uuc,ffc,nnc) = layerer(uu,cc,f)
        decomper_log.info(">>> slicing")
        stdout.flush()
        t3 = timer()
        (ac,kkc) = maxd(nnc)
        ffc1 = fudEmpty()
        if ac > repaRounding:
             ffc1 = depends(ffc,kkc)
        decomper_log.info("dependent fud cardinality : %d" % len(ffqq(ffc1)))
        ccc = apply(uuc,ffc1,cc)
        cc1 = trim(reduce(uuc,fder(ffc1),ccc))
        decomper_log.info("derived cardinality : %d" % acard(cc1))
        qq[(ss,ffc1)] = (ccc,cc1)
        zzc = pathsTree(treesPaths(zz) + [nn+[(ss,ffc1)]])
        (mm,cc,ffc,nnc,kkc) = (None,None,None,None,None)
        t4 = timer()
        decomper_log.info("<<< slicing %.3fs" % (t4-t3))
        stdout.flush()
        return decomp(uuc,zzc,qq,f+1)
    if wmax < 0 or lmax < 0 or xmax < 0 or omax < 0 or bmax < 0 or mmax < 1 or umax < 0 or pmax < 0:
        return None
    if size(aa) == 0 or mult < 1:
        return None
    if not (vars(aa).issubset(uvars(uu)) and vv.issubset(vars(aa))):
        return None
    decomper_log.info(">>> decomper")
    t1 = timer()
    x1 = decomp(uu,emptyTree(),sdict(),1)
    decomper_log.info("nodes: %d" % len(treesNodes(dfzz(x1[1]))))
    t2 = timer()
    decomper_log.info("<<< decomper repa %.3fs" % (t2 - t1))
    stdout.flush()
    return x1

