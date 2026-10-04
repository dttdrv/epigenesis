/* SPDX-License-Identifier: Apache-2.0
 * A local developmental program using the separately licensed CX3D substrate.
 */

import ini.cx3d.Param;
import ini.cx3d.cells.Cell;
import ini.cx3d.cells.CellFactory;
import ini.cx3d.cells.CellModule;
import ini.cx3d.localBiology.AbstractLocalBiologyModule;
import ini.cx3d.localBiology.CellElement;
import ini.cx3d.localBiology.NeuriteElement;
import ini.cx3d.physics.PhysicalCylinder;
import ini.cx3d.physics.PhysicalObject;
import ini.cx3d.physics.PhysicalSphere;
import ini.cx3d.physics.Substance;
import ini.cx3d.simulations.ECM;
import ini.cx3d.simulations.Scheduler;
import ini.cx3d.synapses.Excrescence;

import java.awt.Color;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.Collections;
import java.util.IdentityHashMap;
import java.util.Map;
import java.util.Properties;
import java.util.Set;

import static ini.cx3d.utilities.Matrix.*;


public final class Construction {
    private static final ECM world = ECM.getInstance();
    private static final Properties inputs = new Properties();
    private static BufferedWriter trace;
    private static int step;
    private static int records;
    private static int synapses;
    private static int protrusions;
    private static int retainedSitesOutsideSegments;
    private static java.util.Random contactRandom;
    private static final int MAX_RECORDS = 200_000;
    private static final Set<NeuriteElement> decorated = Collections.newSetFromMap(new IdentityHashMap<>());
    private static final Map<Excrescence,Integer> linked = new IdentityHashMap<>();

    private static double value(String name) {
        double result = Double.parseDouble(inputs.getProperty(name));
        if (!Double.isFinite(result)) throw new IllegalArgumentException("nonfinite input: " + name);
        return result;
    }

    private static String vector(double[] values) {
        StringBuilder result = new StringBuilder("[");
        for (int i = 0; i < values.length; i++) {
            if (!Double.isFinite(values[i])) throw new IllegalStateException("nonfinite physical state");
            if (i != 0) result.append(',');
            result.append(Double.toString(values[i]));
        }
        return result.append(']').toString();
    }

    private static void event(String kind, String fields) {
        if (++records > MAX_RECORDS) throw new IllegalStateException("construction trace budget exhausted");
        try {
            trace.write("{\"event\":\"" + kind + "\",\"step\":" + step + ",\"time\":" + world.getECMtime() + "," + fields + "}\n");
        } catch (IOException failure) {
            throw new java.io.UncheckedIOException(failure);
        }
    }

    private static void budget() {
        if (world.cellList.size() > value("max_cells") || world.neuriteElementList.size() > value("max_segments"))
            throw new IllegalStateException("construction resource budget exhausted");
    }

    private static final class Progenitor implements CellModule {
        private Cell cell;
        private double factor;
        private double bornAt;

        Progenitor(double inheritedFactor) {
            factor = inheritedFactor;
            bornAt = world.getECMtime();
        }

        public Cell getCell() { return cell; }
        public void setCell(Cell next) { cell = next; }
        public boolean isCopiedWhenCellDivides() { return true; }
        public Progenitor getCopy() { return new Progenitor(factor); }

        public void run() {
            budget();
            if (world.getECMtime() <= bornAt) return;
            PhysicalSphere soma = cell.getSomaElement().getPhysicalSphere();
            if (factor > value("division_threshold")) {
                if (soma.getVolume() >= value("division_volume")) {
                    int parent = cell.getID();
                    double before = factor;
                    double volume = soma.getVolume();
                    factor *= 0.5;
                    bornAt = world.getECMtime();
                    Cell daughter = cell.divide(1.0);
                    event("division", "\"parent\":" + parent + ",\"children\":[" + cell.getID() + "," + daughter.getID() + "],"
                        + "\"factor_before\":" + before + ",\"factor_after\":" + factor + ",\"volume_before\":" + volume
                        + ",\"volumes\":" + vector(new double[]{soma.getVolume(),daughter.getSomaElement().getPhysicalSphere().getVolume()}));
                    budget();
                } else {
                    soma.changeVolume(value("volume_rate"));
                }
                return;
            }
            cell.removeCellModule(this);
            NeuriteElement axon = cell.getSomaElement().extendNewNeurite(normalize(randomNoise(1.0,3)));
            axon.setIsAnAxon(true);
            axon.getPhysicalCylinder().setDiameter(value("diameter"));
            axon.addLocalBiologyModule(new Growth(true));
            NeuriteElement dendrite = cell.getSomaElement().extendNewNeurite(normalize(randomNoise(1.0,3)));
            dendrite.setIsAnAxon(false);
            dendrite.getPhysicalCylinder().setDiameter(value("diameter"));
            dendrite.addLocalBiologyModule(new Growth(false));
            event("differentiate", "\"cell\":" + cell.getID() + ",\"factor\":" + factor
                + ",\"soma\":" + soma.getID()
                + ",\"axon\":" + axon.getPhysicalCylinder().getID() + ",\"dendrite\":" + dendrite.getPhysicalCylinder().getID());
        }
    }

    public static double[] growthDirection(double[] direction, double[] gradient, double[] noise, double gain) {
        double[] cue = norm(gradient) == 0.0 ? new double[]{0,0,0} : normalize(gradient);
        double[] combined = add(direction, scalarMult(gain,cue), noise);
        return norm(combined) == 0.0 ? direction.clone() : normalize(combined);
    }

    private static final class Growth extends AbstractLocalBiologyModule {
        private final boolean axon;
        private final double bornAt;
        private double[] direction;

        Growth(boolean isAxon) { axon = isAxon; bornAt = world.getECMtime(); }
        public Growth getCopy() { return new Growth(axon); }
        public boolean isCopiedWhenNeuriteBranches() { return true; }
        public boolean isDeletedAfterNeuriteHasBifurcated() { return true; }
        public void setCellElement(CellElement element) {
            super.setCellElement(element);
            if (!element.isANeuriteElement()) throw new IllegalArgumentException("growth needs a neurite");
            direction = normalize(element.getPhysical().getAxis());
        }

        public void run() {
            budget();
            if (world.getECMtime() <= bornAt) return;
            NeuriteElement process = (NeuriteElement) cellElement;
            PhysicalCylinder cylinder = process.getPhysicalCylinder();
            double diameter = cylinder.getDiameter();
            if (diameter <= value("stop_diameter")) {
                process.removeLocalBiologyModule(this);
                event("stop", "\"segment\":" + cylinder.getID() + ",\"diameter\":" + diameter);
                return;
            }
            double[] gradient = cylinder.getExtracellularGradient("guidance");
            double[] priorDirection = direction;
            double[] noise = randomNoise(value("noise"),3);
            direction = growthDirection(direction,gradient,noise,value("guidance_gain"));
            double speed = value(axon ? "axon_speed" : "dendrite_speed");
            event("growth", "\"cell\":" + process.getCell().getID() + ",\"segment\":" + cylinder.getID()
                + ",\"position\":" + vector(cylinder.getMassLocation()) + ",\"direction\":" + vector(direction)
                + ",\"prior_direction\":" + vector(priorDirection) + ",\"noise\":" + vector(noise) + ",\"axon\":" + axon
                + ",\"gradient\":" + vector(gradient) + ",\"speed\":" + speed + ",\"dt\":" + Param.SIMULATION_TIME_STEP);
            cylinder.movePointMass(speed,direction);
            cylinder.setDiameter(Math.max(value("stop_diameter"),diameter-value("taper")*Param.SIMULATION_TIME_STEP));
            double probability = -Math.expm1(-value("branch_rate")*Param.SIMULATION_TIME_STEP);
            if (ECM.getRandomDouble() < probability && cylinder.getDiameter()/Math.sqrt(2) > value("stop_diameter")) {
                double branchDiameter = cylinder.getDiameter()/Math.sqrt(2);
                NeuriteElement[] children = process.bifurcate();
                for (NeuriteElement child : children) {
                    child.setIsAnAxon(axon);
                    child.getPhysicalCylinder().setDiameter(branchDiameter);
                }
                event("branch", "\"cell\":" + process.getCell().getID() + ",\"parent\":" + cylinder.getID()
                    + ",\"children\":[" + children[0].getPhysicalCylinder().getID() + "," + children[1].getPhysicalCylinder().getID()
                    + "],\"diameter_before\":" + cylinder.getDiameter() + ",\"diameter_after\":" + branchDiameter);
                budget();
            }
        }
    }

    public static int protrusionCount(double length, double spacing, int remaining) {
        double count = length/spacing;
        if (!Double.isFinite(count) || length < 0 || spacing <= 0 || remaining < 0 || Math.round(count) > remaining)
            throw new IllegalStateException("construction protrusion budget exhausted");
        return (int)Math.round(count);
    }

    private static void connect() {
        if (value("recognition") == 0.0) return;
        for (NeuriteElement process : world.neuriteElementList) {
            PhysicalCylinder cylinder = process.getPhysicalCylinder();
            if (!cylinder.isTerminal() && decorated.add(process)) {
                int count = protrusionCount(cylinder.getActualLength(),value("contact_spacing"),(int)value("max_protrusions")-protrusions);
                int before = cylinder.getExcrescences().size();
                if (process.isAnAxon()) process.makeBoutons(value("contact_spacing"));
                else process.makeSpines(value("contact_spacing"));
                if (cylinder.getExcrescences().size()-before != count) throw new IllegalStateException("protrusion count mismatch");
                protrusions += count;
                for (Excrescence protrusion : cylinder.getExcrescences()) protrusion.setLength(value("contact_reach")/2);
            }
        }
        for (NeuriteElement process : world.neuriteElementList) {
            if (process.isAnAxon()) process.synapseBetweenExistingBS(1.0);
        }
        for (NeuriteElement process : world.neuriteElementList) {
            if (!process.isAnAxon()) continue;
            for (Excrescence bouton : process.getPhysicalCylinder().getExcrescences()) {
                Excrescence spine = bouton.getEx();
                if (spine == null || linked.containsKey(bouton)) continue;
                if (spine.getEx() != bouton || bouton.getType() != Excrescence.BOUTON || spine.getType() != Excrescence.SPINE)
                    throw new IllegalStateException("invalid physical synapse");
                int pre = process.getCell().getID();
                int post = spine.getPo().getCellElement().getCell().getID();
                if (pre == post) throw new IllegalStateException("unexpected autapse");
                synapses++;
                linked.put(bouton,synapses);
                event("synapse", "\"contact\":" + synapses + ",\"pre\":" + pre + ",\"post\":" + post + ",\"pre_segment\":" + bouton.getPo().getID()
                    + ",\"post_segment\":" + spine.getPo().getID() + ",\"pre_position\":" + vector(bouton.getProximalEnd())
                    + ",\"post_position\":" + vector(spine.getProximalEnd()) + ",\"reach\":" + (bouton.getLength()+spine.getLength()));
            }
        }
    }

    private static void connections() {
        int count = 0;
        for (NeuriteElement process : world.neuriteElementList) {
            if (!process.isAnAxon()) continue;
            for (Excrescence bouton : process.getPhysicalCylinder().getExcrescences()) {
                if (!linked.containsKey(bouton)) continue;
                Excrescence spine = bouton.getEx();
                if (spine == null || spine.getEx() != bouton) throw new IllegalStateException("lost reciprocal contact");
                double preAxial = bouton.getPositionOnPO()[0];
                double postAxial = spine.getPositionOnPO()[0];
                double preLength = ((PhysicalCylinder)bouton.getPo()).getActualLength();
                double postLength = ((PhysicalCylinder)spine.getPo()).getActualLength();
                if (preAxial < 0 || preAxial > preLength) retainedSitesOutsideSegments++;
                if (postAxial < 0 || postAxial > postLength) retainedSitesOutsideSegments++;
                event("connection", "\"contact\":" + linked.get(bouton) + ",\"pre\":" + process.getCell().getID()
                    + ",\"post\":" + spine.getPo().getCellElement().getCell().getID()
                    + ",\"pre_segment\":" + bouton.getPo().getID() + ",\"post_segment\":" + spine.getPo().getID()
                    + ",\"pre_axial\":" + preAxial + ",\"post_axial\":" + postAxial
                    + ",\"pre_length\":" + preLength + ",\"post_length\":" + postLength);
                count++;
            }
        }
        if (count != synapses) throw new IllegalStateException("contact retention disagrees with formation history");
    }

    private static void snapshot() {
        for (Cell cell : world.cellList) {
            PhysicalSphere soma = cell.getSomaElement().getPhysicalSphere();
            event("cell", "\"cell\":" + cell.getID() + ",\"position\":" + vector(soma.getMassLocation())
                + ",\"soma\":" + soma.getID() + ",\"volume\":" + soma.getVolume());
        }
        for (PhysicalCylinder cylinder : world.physicalCylinderList) {
            PhysicalObject parent = cylinder.getMother();
            event("segment", "\"segment\":" + cylinder.getID() + ",\"cell\":" + cylinder.getCellElement().getCell().getID()
                + ",\"parent\":" + parent.getID() + ",\"parent_is_soma\":" + (parent instanceof PhysicalSphere)
                + ",\"position\":" + vector(cylinder.getMassLocation()) + ",\"proximal\":" + vector(subtract(cylinder.getMassLocation(),cylinder.getSpringAxis()))
                + ",\"diameter\":" + cylinder.getDiameter() + ",\"axon\":" + cylinder.getNeuriteElement().isAnAxon());
        }
    }

    public static void main(String[] args) throws Exception {
        if (args.length != 3) throw new IllegalArgumentException("expected input, trace, and summary paths");
        try (java.io.Reader input = Files.newBufferedReader(Path.of(args[0]))) { inputs.load(input); }
        Param.SIMULATION_TIME_STEP = value("dt");
        ECM.setRandomSeed(Long.parseLong(inputs.getProperty("seed")));
        contactRandom = new java.util.Random(new java.util.SplittableRandom(Long.parseLong(inputs.getProperty("seed"))).split().nextLong());
        double[] origin = {value("origin_x"),value("origin_y"),value("origin_z")};
        for (int i = 0; i < 18; i++) world.getPhysicalNodeInstance(add(origin,randomNoise(1000,3)));
        world.addArtificialGaussianConcentrationZ(new Substance("guidance",Color.RED),1.0,origin[2]+value("cue"),40.0);
        Cell founder = CellFactory.getCellInstance(origin);
        founder.addCellModule(new Progenitor(value("division_signal")));
        Scheduler.setPrintCurrentECMTime(false);
        PhysicalCylinder.subdivisionObserver = (distal,proximal) -> event("subdivide", "\"distal\":" + distal.getID()
            + ",\"proximal_segment\":" + proximal.getID() + ",\"parent\":" + proximal.getMother().getID());
        PhysicalCylinder.mergeObserver = (distal,proximal) -> event("merge", "\"distal\":" + distal.getID()
            + ",\"proximal_segment\":" + proximal.getID() + ",\"parent\":" + distal.getMother().getID());
        world.canRun.release();
        try (BufferedWriter writer = Files.newBufferedWriter(Path.of(args[1]),StandardOpenOption.CREATE_NEW)) {
            trace = writer;
            event("founder", "\"cell\":" + founder.getID() + ",\"position\":" + vector(origin)
                + ",\"factor\":" + value("division_signal"));
            snapshot();
            for (step = 1; step <= (int)value("steps"); step++) {
                if (step == (int)value("cue_switch_step")) {
                    world.addArtificialGaussianConcentrationZ("guidance",1.0,origin[2]+value("cue_after"),40.0);
                    event("environment", "\"cue\":" + value("cue_after"));
                }
                Scheduler.simulateOneStep();
                budget();
                ECM.withRandom(contactRandom, Construction::connect);
                snapshot();
            }
            step--;
            connections();
            double length = 0;
            for (PhysicalCylinder cylinder : world.physicalCylinderList) length += cylinder.getActualLength();
            String summary = "{\"initial_cells\":1,\"cells\":" + world.cellList.size() + ",\"segments\":" + world.neuriteElementList.size()
                + ",\"synapses\":" + synapses + ",\"length\":" + length + ",\"time\":" + world.getECMtime()
                + ",\"retained_sites_outside_segments\":" + retainedSitesOutsideSegments
                + ",\"physics\":" + Scheduler.runPhyics + ",\"diffusion\":" + Scheduler.runDiffusion + "}\n";
            Files.writeString(Path.of(args[2]),summary,StandardOpenOption.CREATE_NEW);
        }
    }
}
