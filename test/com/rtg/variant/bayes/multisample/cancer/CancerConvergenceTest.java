/*
 * Copyright (c) 2014. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */

package com.rtg.variant.bayes.multisample.cancer;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import com.rtg.relation.GenomeRelationships;
import com.rtg.relation.Relationship.RelationshipType;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.Variant;
import com.rtg.variant.VariantOutputLevel;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantParamsBuilder;
import com.rtg.variant.bayes.EvidenceInterface;
import com.rtg.variant.bayes.Model;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.multisample.ComplexCallerTest;
import com.rtg.variant.bayes.multisample.HaploidDiploidHypotheses;
import com.rtg.variant.bayes.multisample.MultisampleJointCaller;
import com.rtg.variant.bayes.snp.DescriptionSnp;
import com.rtg.variant.bayes.snp.EvidenceQ;
import com.rtg.variant.bayes.snp.HypothesesNone;
import com.rtg.variant.bayes.snp.HypothesesSnp;
import com.rtg.variant.bayes.snp.StatisticsSnp;
import com.rtg.variant.util.VariantUtils;
import com.rtg.variant.util.arithmetic.LogPossibility;

import htsjdk.samtools.SAMFileHeader;
import junit.framework.TestCase;

/**
 */
public class CancerConvergenceTest extends TestCase {

  private SomaticCallerConfiguration getConfig(final double contamination) throws IOException, InvalidParamsException {
    final GenomeRelationships genomeRelationships = new GenomeRelationships();
    genomeRelationships.addGenome("cancer");
    genomeRelationships.addGenome("TEST");
    genomeRelationships.addRelationship(RelationshipType.ORIGINAL_DERIVED, "TEST", "cancer").setProperty("contamination", String.valueOf(contamination));

    final SAMFileHeader uber = ComplexCallerTest.makeHeaderWithSamples("TEST", "cancer");
    final VariantParamsBuilder b = new VariantParamsBuilder();
    b.genomePriors(GenomePriorParams.builder().create());
    b.genomeRelationships(genomeRelationships);
    b.machineErrorName("illumina");
    b.lohPrior(1e-20);
    b.uberHeader(uber);
    final VariantParams p = b.create();
    return new SomaticCallerConfiguration.Configurator().getConfig(p, null);
  }

  private MultisampleJointCaller getCaller(final SomaticCallerConfiguration jointConfig) {
    return jointConfig.getJointCaller();
  }

  private static final int MAPQ = 20;
  private static final double MAPQ_PROB = VariantUtils.phredToProb(MAPQ);
  private static final int PHRED = 20;
  private static final double PHREAD_PROB = VariantUtils.phredToProb(PHRED);

  private HypothesesSnp getHypotheses(final int refNt) {
    return new HypothesesSnp(LogPossibility.SINGLETON, GenomePriorParams.builder().create(), false, refNt);
  }
  private HypothesesSnp getHypothesesHaploid(final int refNt) {
    return new HypothesesSnp(LogPossibility.SINGLETON, GenomePriorParams.builder().create(), true, refNt);
  }

  private int randomNot(final Random r, final int v) {
    while (true) {
      final int q = r.nextInt(4);
      if (q != v) {
        return q;
      }
    }
  }

  public int[] checkHeterozygousConvergenceWithNoise(final double noiseLevel, final double contamination) throws IOException, InvalidParamsException {
    Diagnostic.setLogStream();
    final Random r = new Random(3231);
    final int maxCoverage = 100;
    final int[] correct = new int[maxCoverage];
    final SomaticCallerConfiguration jointConfig = getConfig(contamination);
    final MultisampleJointCaller jointCaller = getCaller(jointConfig);
    final GenomePriorParams gp = GenomePriorParams.builder().create();
    final VariantParams params = VariantParams.builder().callLevel(VariantOutputLevel.ALL).genomePriors(gp).create();
    final ModelCancerFactory mf = new ModelCancerFactory(params.genomePriors(), contamination, false);

    for (int j = 0; j < 100; j++) {
      final byte[] ref = {(byte) (1 + r.nextInt(4))}; // random reference A, C, G, T
      final HypothesesSnp hypotheses = getHypotheses(ref[0] - 1);
      final HypothesesSnp hypothesesHaploid = getHypothesesHaploid(ref[0] - 1);
      final ModelInterface<?> normalModel = new Model<>(hypotheses, new StatisticsSnp(hypotheses.description()));
      final ModelInterface<?> cancerModel = mf.make(ref[0] - 1);

      // Choose random alleles for normal and cancer, such that there is always an explanation.
      int normal1, normal2, cancer1, cancer2;
      do {
        normal1 = r.nextInt(4);
        normal2 = r.nextInt(4);
        cancer1 = r.nextInt(4);
        cancer2 = r.nextInt(4);
      } while ((normal1 == cancer1 && normal2 == cancer2) || (normal1 == cancer2 && normal2 == cancer1) || (cancer1 == cancer2 && cancer1 == ref[0] - 1));


      //System.out.println("normal = " + normal1 + ":" + normal2 + " cancer = " + cancer1 + ":" + cancer2);

      // Random generation below is not perfect in the simulation of noise. It tries to make different nucleotides
      // based on selecting an allele and then changing that particular allele.

      for (int k = 0; k < maxCoverage; k++) {
        final boolean noise = r.nextDouble() < noiseLevel;
        final EvidenceInterface evidence = new EvidenceQ(DescriptionSnp.SINGLETON, noise ? randomNot(r, r.nextBoolean() ? normal1 : normal2) : r.nextBoolean() ? normal1 : normal2, 0, 0, MAPQ_PROB, PHREAD_PROB, true, false, false, false);
        normalModel.increment(evidence);
        if (r.nextDouble() < contamination) {
          final EvidenceInterface evidence1 = new EvidenceQ(DescriptionSnp.SINGLETON, noise ? randomNot(r, r.nextBoolean() ? normal1 : normal2) : r.nextBoolean() ? normal1 : normal2, 0, 0, MAPQ_PROB, PHREAD_PROB, true, false, false, false);
          cancerModel.increment(evidence1);
        } else {
          final EvidenceInterface evidence1 = new EvidenceQ(DescriptionSnp.SINGLETON, noise ? randomNot(r, r.nextBoolean() ? cancer1 : cancer2) : r.nextBoolean() ? cancer1 : cancer2, 0, 0, MAPQ_PROB, PHREAD_PROB, true, false, false, false);
          cancerModel.increment(evidence1);
        }
        final List<ModelInterface<?>> models = new ArrayList<>();
        models.add(normalModel);
        models.add(cancerModel);
        final Variant v = jointCaller.makeCall("refName", 0, 1, ref, models, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, hypothesesHaploid, hypotheses));
        if (v != null) {
          final String normal = v.getSample(0).getName();
          final String cancer = v.getSample(1).getName();
          if (!"NONE".equals(normal) && !"NONE".equals(cancer)) {
            //System.out.println("normal=" + normal + " cancer=" + cancer);
            final int n1 = "ACGT".indexOf(normal.charAt(0));
            final int n2 = normal.contains(":") ? "ACGT".indexOf(normal.charAt(2)) : n1;
            final int c1 = "ACGT".indexOf(cancer.charAt(0));
            final int c2 = cancer.contains(":") ? "ACGT".indexOf(cancer.charAt(2)) : c1;
            if (((n1 == normal1 && n2 == normal2) || (n1 == normal2 && n2 == normal1)) && ((c1 == cancer1 && c2 == cancer2) || (c1 == cancer2 && c2 == cancer1))) {
              correct[k]++;
            }
          }
        }
      }
    }
    return correct;
  }

  public void testHeterozygousConvergenceWithNoiseNoContam() throws IOException, InvalidParamsException {
    final int[] correct00 = checkHeterozygousConvergenceWithNoise(0, 0);
    final int[] correct01 = checkHeterozygousConvergenceWithNoise(0.1, 0);
    final int[] correct02 = checkHeterozygousConvergenceWithNoise(0.2, 0);

    //    for (int k = 0; k < correct00.length; k++) {
    //      System.out.println(k + " " + correct00[k] + " " + correct01[k] + " " + correct02[k]);
    //    }
    // We have a curve showing how often we find the answer.  The following is essentially
    // a very crude hull under the curve.
    for (int k = 7; k < correct00.length; k++) {
      assertTrue(correct00[k] >= k);
      assertTrue(correct01[k] >= k);
      assertTrue(correct02[k] >= k - 3);
    }
  }

  public void testHeterozygousConvergenceWithNoiseContam() throws IOException, InvalidParamsException {
    final int[] correct00 = checkHeterozygousConvergenceWithNoise(0, 0.3);
    final int[] correct01 = checkHeterozygousConvergenceWithNoise(0.1, 0.3);
    final int[] correct02 = checkHeterozygousConvergenceWithNoise(0.2, 0.3);

    //    for (int k = 0; k < correct00.length; k++) {
    //      System.out.println(k + " " + correct00[k] + " " + correct01[k] + " " + correct02[k]);
    //    }
    // We have a curve showing how often we find the answer.  The following is essentially
    // a very crude hull under the curve.
    for (int k = 7; k < correct00.length; k++) {
      assertTrue(correct00[k] >= k);
      assertTrue(correct01[k] >= k);
      assertTrue(correct02[k] >= k - 20);
    }
  }

  // Uncomment the following to produce table of numbers showing convergence varies with contamination
  //  public void testHeterozygousConvergenceWithNoiseContamGraph(), InvalidParamsException {
  //    System.out.println("plot \"convergence\" using 1:2 with linesp title \"0%\", \"convergence\" using 1:3 with linesp title \"25%\", \"convergence\" using 1:4 with linesp title \"50%\", \"convergence\" using 1:5 with linesp title \"75%\",  \"convergence\" using 1:6 with linesp title \"90%\"");
  //    final double error = 0.1;
  //    final ArrayList<int[]> res = new ArrayList<>();
  //    for (final double contamination : new double[] {0, 0.25, 0.5, 0.75, 0.9}) {
  //      res.add(checkHeterozygousConvergenceWithNoise(error, contamination));
  //    }
  //    for (int k = 0; k < res.get(0).length; k++) {
  //      final StringBuilder sb = new StringBuilder();
  //      sb.append(k);
  //      for (int j = 0; j < res.size(); j++) {
  //        sb.append(" ").append(res.get(j)[k]);
  //      }
  //      System.out.println(sb.toString());
  //    }
  //  }
}
