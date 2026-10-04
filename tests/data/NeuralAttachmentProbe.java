/* SPDX-License-Identifier: Apache-2.0 */
import ini.cx3d.Param;
import ini.cx3d.cells.CellFactory;
import ini.cx3d.physics.*;
import ini.cx3d.simulations.ECM;
import ini.cx3d.synapses.*;
import java.util.*;
import static ini.cx3d.utilities.Matrix.*;

public class NeuralAttachmentProbe {
    static final List<Excrescence> sites = new ArrayList<>();
    static final List<PhysicalBond> bonds = new ArrayList<>();
    static final Map<PhysicalBond,double[]> positions = new IdentityHashMap<>();
    static final Map<PhysicalBond,double[]> coordinates = new IdentityHashMap<>();
    static final Map<PhysicalBond,Double> rest = new IdentityHashMap<>();

    static PhysicalCylinder cylinder(double y) {
        return CellFactory.getCellInstance(new double[]{0,y,0}).getSomaElement()
            .extendNewNeurite(new double[]{1,0,0}).getPhysicalCylinder();
    }

    static double[] endpoint(PhysicalBond bond, PhysicalObject other) {
        return bond.getFirstPhysicalObject()==other ? bond.getSecondEndLocation() : bond.getFirstEndLocation();
    }

    static void attach(PhysicalCylinder owner, PhysicalCylinder other, double fraction) {
        double[] point = {fraction*owner.getActualLength(),0};
        PhysicalBouton bouton = new PhysicalBouton(owner,point,2);
        PhysicalSpine spine = new PhysicalSpine(other,new double[]{.5,Math.PI},2);
        owner.addExcrescence(bouton); other.addExcrescence(spine);
        assert bouton.synapseWith(spine,true);
        sites.add(bouton);
        // includes another bond sharing this exact coordinate array.
        new PhysicalBond(owner,point,other,new double[]{.5,0},1,1);
    }

    static void deformation(String mode) {
        PhysicalCylinder owner = cylinder(0), other = cylinder(100);
        owner.movePointMass(9/Param.SIMULATION_TIME_STEP,new double[]{1,0,0});
        List<PhysicalCylinder> affected = new ArrayList<>(); affected.add(owner);
        if (mode.equals("passive-branch")) {
            owner.getNeuriteElement().bifurcate();
            affected.add(owner.getDaughterLeft()); affected.add(owner.getDaughterRight());
        }
        Map<PhysicalCylinder,Double> lengths = new IdentityHashMap<>();
        Map<PhysicalCylinder,Double> volumes = new IdentityHashMap<>();
        Map<double[],Double> axial = new IdentityHashMap<>();
        Map<PhysicalBond,Double> rests = new IdentityHashMap<>();
        for (PhysicalCylinder cylinder : affected) {
            double length = cylinder.getActualLength();
            lengths.put(cylinder,length); volumes.put(cylinder,cylinder.getVolume());
            IntracellularSubstance substance = new IntracellularSubstance("strain",0,0);
            substance.setQuantity(10); cylinder.getIntracellularSubstances().put("strain",substance);
            for (double fraction : new double[]{0,.5,.999,1}) {
                double[] point = {fraction*length,.3};
                PhysicalBouton site = new PhysicalBouton(cylinder,point,2);
                cylinder.addExcrescence(site); axial.put(point,point[0]);
                for (int i=0;i<2;i++) {
                    double[] remote = {.5,0};
                    double rest = distance(cylinder.transformCoordinatesPolarToGlobal(point),other.transformCoordinatesPolarToGlobal(remote));
                    new PhysicalBond(cylinder,point,other,remote,rest,1);
                }
            }
            // independent endpoints, including the distal point mass, must also move.
            new PhysicalBond(cylinder,other);
            new PhysicalBond(other,other.transformCoordinatesGlobalToPolar(other.getMassLocation()),
                cylinder,new double[]{.5*length,0},100,0);
            for (PhysicalBond bond : cylinder.getPhysicalBonds()) {
                double[] point = bond.getPositionOnObjectInLocalCoord(cylinder);
                axial.put(point,point[0]); rests.put(bond,bond.getRestingLength());
            }
        }
        if (mode.equals("passive-invalid-length")) {
            try {
                java.lang.reflect.Method update = PhysicalCylinder.class.getDeclaredMethod("updateAfterPassiveDeformation");
                update.setAccessible(true);
                for (double invalid : new double[]{0,-1,Double.NaN,Double.POSITIVE_INFINITY}) {
                    owner.setActualLength(invalid);
                    try { update.invoke(owner); throw new AssertionError("invalid prior length accepted"); }
                    catch (java.lang.reflect.InvocationTargetException failure) {
                        assert failure.getCause() instanceof IllegalStateException;
                    }
                    for (Excrescence site : owner.getExcrescences())
                        assert site.getPositionOnPO()[0]==axial.get(site.getPositionOnPO());
                    owner.setActualLength(lengths.get(owner));
                    update.invoke(owner);
                }
            } catch (ReflectiveOperationException failure) { throw new AssertionError(failure); }
            return;
        }
        boolean growth = mode.equals("active-growth");
        owner.setRestingLengthForDesiredTension(mode.equals("passive-stretch")?-1:10);
        if (growth) owner.movePointMass(1/Param.SIMULATION_TIME_STEP,new double[]{1,0,0});
        else if (mode.equals("passive-soma")) ((PhysicalSphere)owner.getMother()).runPhysics();
        else owner.runPhysics();
        int changed = 0;
        for (PhysicalCylinder cylinder : affected) {
            double before = lengths.get(cylinder), after = cylinder.getActualLength();
            if (Math.abs(before-after)>1e-10) changed++;
            for (Excrescence site : cylinder.getExcrescences()) {
                double[] point = site.getPositionOnPO();
                double expected = growth?axial.get(point):(axial.get(point)/before)*after;
                assert Math.abs(point[0]-expected)<1e-12 : "attachment did not follow declared material motion";
                assert point[0]>=0 && point[0]<=after : "passive strain moved valid site outside shaft";
                assert point[1]==.3 && site.getPo()==cylinder : "attachment angle or owner changed";
                long aliases = cylinder.getPhysicalBonds().stream()
                    .filter(b->b.getPositionOnObjectInLocalCoord(cylinder)==point).count();
                assert aliases==2 : "shared coordinate identity changed";
            }
            for (PhysicalBond bond : cylinder.getPhysicalBonds()) {
                double[] point = bond.getPositionOnObjectInLocalCoord(cylinder);
                double expected = growth?axial.get(point):(axial.get(point)/before)*after;
                assert Math.abs(point[0]-expected)<1e-12 : "independent bond coordinate did not follow strain";
                assert bond.getRestingLength()==rests.get(bond) : "strain changed bond rest length";
                assert Collections.frequency(cylinder.getPhysicalBonds(),bond)==1;
                assert Collections.frequency(other.getPhysicalBonds(),bond)==1;
                if (!growth && axial.get(point)==before && point.length==3 && point[2]==0)
                    assert distance(endpoint(bond,other),cylinder.getMassLocation())<1e-12;
            }
            assert cylinder.getIntracellularSubstance("strain").getQuantity()==10;
            assert Math.abs(cylinder.getVolume()/volumes.get(cylinder)-after/before)<1e-12;
        }
        assert changed==affected.size() : "probe did not deform every intended segment";
        if (mode.equals("passive-stretch")) assert owner.getActualLength()>lengths.get(owner);
    }

    static void invalidDeformation(String mode) {
        PhysicalCylinder owner = cylinder(0);
        owner.movePointMass(9/Param.SIMULATION_TIME_STEP,new double[]{1,0,0});
        List<PhysicalCylinder> affected = new ArrayList<>(); affected.add(owner);
        if (mode.equals("passive-collapsed-length")) {
            owner.addExcrescence(new PhysicalBouton(owner,new double[]{0,0},2));
            owner.setMassLocation(owner.proximalEnd());
            owner.setActualLength(0);
            try {
                java.lang.reflect.Method update = PhysicalCylinder.class.getDeclaredMethod("updateAfterPassiveDeformation");
                update.setAccessible(true);
                try { update.invoke(owner); throw new AssertionError("collapsed material length accepted"); }
                catch (java.lang.reflect.InvocationTargetException failure) {
                    assert failure.getCause() instanceof IllegalStateException;
                }
            } catch (ReflectiveOperationException failure) { throw new AssertionError(failure); }
        } else {
            PhysicalCylinder[] children = owner.bifurcateCylinder(.125,new double[]{1,0,0},new double[]{1,1,0});
            affected.addAll(Arrays.asList(children));
            children[0].addExcrescence(new PhysicalBouton(children[0],new double[]{.1,0},2));
            Param.SIMULATION_TIME_STEP = .125;
            double[] next = owner.getMassLocation().clone(), axis = owner.getSpringAxis();
            for (int i=0;i<3;i++) next[i] += axis[i]/owner.getActualLength()*Param.SIMULATION_TIME_STEP;
            children[0].setMassLocation(next);
            try {
                java.lang.reflect.Method update = PhysicalCylinder.class.getDeclaredMethod("updateDependentPhysicalVariables");
                update.setAccessible(true); update.invoke(children[0]);
            } catch (ReflectiveOperationException failure) { throw new AssertionError(failure); }
            assert children[0].getActualLength()>=.1 : "probe starts with invalid attachment geometry";
            children[0].setRestingLengthForDesiredTension(0);
            children[1].setRestingLengthForDesiredTension(0);
            owner.setRestingLengthForDesiredTension(-1);
            try { owner.runPhysics(); throw new AssertionError("invalid child length accepted"); }
            catch (IllegalStateException expected) { }
        }
        for (PhysicalCylinder cylinder : affected) {
            java.util.concurrent.locks.ReentrantReadWriteLock lock =
                (java.util.concurrent.locks.ReentrantReadWriteLock)cylinder.getRwLock();
            assert lock.getReadHoldCount()==0 && lock.getWriteHoldCount()==0 : "deformation failure leaked a lock";
        }
    }

    public static void main(String[] args) {
        ECM world = ECM.getInstance(); ECM.setRandomSeed(7);
        for (int i=0;i<18;i++) world.getPhysicalNodeInstance(randomNoise(1000,3));
        if (args[0].equals("passive-collapsed-length") || args[0].equals("passive-child-failure")) {
            invalidDeformation(args[0]); return;
        }
        if (args[0].startsWith("passive-") || args[0].equals("active-growth")) {
            deformation(args[0]); return;
        }
        PhysicalCylinder distal = cylinder(0), other = cylinder(3), proximal = null;
        boolean split = args[0].equals("split");
        boolean retain = args[0].equals("tapered") || args[0].equals("rotated");
        distal.movePointMass((split?19:5)/Param.SIMULATION_TIME_STEP,new double[]{1,0,0});
        IntracellularSubstance substance = new IntracellularSubstance("remesh",1,0);
        substance.setQuantity(10); distal.getIntracellularSubstances().put("remesh",substance);
        if (!split) {
            double maximum = Param.NEURITE_MAX_LENGTH;
            Param.NEURITE_MAX_LENGTH = 5;
            distal.runDiscretization();
            Param.NEURITE_MAX_LENGTH = maximum;
            distal.bifurcateCylinder(1,new double[]{1,1,0},new double[]{1,-1,0});
            proximal = (PhysicalCylinder)distal.getMother();
        }
        if (args[0].equals("tapered")) proximal.setDiameter(2*distal.getDiameter());
        if (args[0].equals("rotated")) {
            double[] y = proximal.getYAxis();
            proximal.setYAxis(proximal.getZAxis());
            proximal.setZAxis(scalarMult(-1,y));
        }
        if (!args[0].equals("proximal-only")) {
            for (double fraction : new double[]{0,.5,.9,.95,1}) attach(distal,other,fraction);
        }
        if (proximal != null && !args[0].equals("distal-only")) attach(proximal,other,.5);
        // independent endpoints without any excrescence must also be remapped.
        new PhysicalBond(distal,new double[]{.95*distal.getActualLength(),0},other,new double[]{.5,0},1,1);
        new PhysicalBond(other,new double[]{.5,0},distal,new double[]{.1*distal.getActualLength(),0},1,1);
        if (proximal != null) new PhysicalBond(proximal,new double[]{.5*proximal.getActualLength(),0},other,new double[]{.5,0},1,1);
        bonds.addAll(other.getPhysicalBonds());
        for (PhysicalBond bond : bonds) {
            positions.put(bond,endpoint(bond,other));
            coordinates.put(bond,bond.getPositionOnObjectInLocalCoord(bond.getOppositePhysicalObject(other)));
            rest.put(bond,bond.getRestingLength());
        }
        double volume = distal.getVolume()+(proximal==null?0:proximal.getVolume());
        double quantity = distal.getIntracellularSubstance("remesh").getQuantity()
            +(proximal==null?0:proximal.getIntracellularSubstance("remesh").getQuantity());
        double length = distal.getActualLength()+(proximal==null?0:proximal.getActualLength());
        distal.runDiscretization();
        List<PhysicalCylinder> owners = new ArrayList<>(); owners.add(distal);
        if (split) owners.add((PhysicalCylinder)distal.getMother());
        if (retain) {
            assert distal.getMother()==proximal && proximal.isStillExisting() : "incompatible geometry was merged";
            owners.add(proximal);
        }
        Set<Excrescence> seen = Collections.newSetFromMap(new IdentityHashMap<>());
        double finalVolume=0, finalQuantity=0, finalLength=0;
        for (PhysicalCylinder owner : owners) {
            finalVolume += owner.getVolume(); finalLength += owner.getActualLength();
            finalQuantity += owner.getIntracellularSubstance("remesh").getQuantity();
            for (Excrescence site : owner.getExcrescences()) {
                assert site.getPo()==owner && seen.add(site) : "wrong or duplicated attachment owner";
            }
        }
        assert seen.size()==sites.size() && seen.containsAll(sites) : "lost attachment identity";
        assert Math.abs(finalVolume-volume)<1e-10 && Math.abs(finalQuantity-quantity)<1e-12
            && Math.abs(finalLength-length)<1e-12 : "remesh changes conserved geometry or substance";
        if (proximal!=null && !retain) assert proximal.getExcrescences().isEmpty() && proximal.getPhysicalBonds().isEmpty();
        for (PhysicalBond bond : bonds) {
            PhysicalObject owner = bond.getOppositePhysicalObject(other);
            assert owners.contains(owner) && owner.isStillExisting() : "bond has stale owner";
            assert Collections.frequency(owner.getPhysicalBonds(),bond)==1;
            assert Collections.frequency(other.getPhysicalBonds(),bond)==1;
            assert bond.getPositionOnObjectInLocalCoord(owner)==coordinates.get(bond) : "coordinate alias lost";
            assert distance(positions.get(bond),endpoint(bond,other))<1e-12 : "straight remesh moved attachment";
            assert bond.getRestingLength()==rest.get(bond) : "remesh changed bond rest length";
            for (Excrescence site : sites) if (site.getPositionOnPO()==coordinates.get(bond)) {
                assert site.getPo()==owner && distance(site.getProximalEnd(),endpoint(bond,other))<1e-12;
            }
        }
    }
}
