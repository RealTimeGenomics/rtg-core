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
package com.rtg.variant.bayes.multisample.family;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 * Test class for all tests in this package.
 *
 */
public class AllTests extends TestSuite {

  public static Test suite() {
    final TestSuite suite = new TestSuite(AllTests.class.getName());
    suite.addTestSuite(AlleleProbabilityHNTest.class);
    suite.addTestSuite(AlleleProbabilityHHTest.class);
    suite.addTestSuite(AlleleSetProbabilityDiploidTest.class);
    suite.addTestSuite(AlleleSetProbabilityHaploidTest.class);
    suite.addTestSuite(BinomialSpecialTest.class);
    suite.addTestSuite(BitSetTest.class);
    suite.addTestSuite(DefaultWeightedLatticeTest.class);
    suite.addTestSuite(DescriptionDiseaseTest.class);
    suite.addTestSuite(DiseasedFamilyCallerTest.class);
    suite.addTestSuite(DiseasedFamilyCallerConfigurationTest.class);
    suite.addTestSuite(DiseasedFamilyPosteriorTest.class);
    suite.addTestSuite(DiseaseTest.class);
    suite.addTestSuite(MendelianAlleleProbabilityNormalTest.class);
    suite.addTestSuite(MendelianAlleleProbabilityDeNovoTest.class);
    suite.addTestSuite(DummySegregationHaploidTest.class);
    suite.addTestSuite(FamilyCliTest.class);
    suite.addTestSuite(FamilyCallerTest.class);
    suite.addTestSuite(FamilyCallerConfigurationTest.class);
    suite.addTestSuite(FamilyNanoTest.class);
    suite.addTestSuite(FamilyPosteriorTest.class);
    suite.addTestSuite(FastDiseasedFamilyPosteriorTest.class);
    suite.addTestSuite(FastFamilyPosteriorTest.class);
    suite.addTestSuite(HypothesesDiseaseTest.class);
    suite.addTestSuite(MendelianAlleleProbabilityCombinerTest.class);
    suite.addTestSuite(MendelianAlleleProbabilityDiploidTest.class);
    suite.addTestSuite(MendelianAlleleProbabilityHDDTest.class);
    suite.addTestSuite(MendelianAlleleProbabilityHDHTest.class);
    suite.addTestSuite(MendelianAlleleProbabilityHHHTest.class);
    suite.addTestSuite(MendelianAlleleProbabilityNHHTest.class);
    suite.addTestSuite(MendelianAlleleProbabilityNHNTest.class);
    suite.addTestSuite(MendelianAlleleProbabilityDiploidDeNovoTest.class);
    suite.addTestSuite(MendelianAlleleProbabilityHDDDeNovoTest.class);
    suite.addTestSuite(MendelianAlleleProbabilityHDHDeNovoTest.class);
    suite.addTestSuite(MendelianAlleleProbabilityHHHDeNovoTest.class);
    suite.addTestSuite(MendelianAlleleProbabilityNHHDeNovoTest.class);
    suite.addTestSuite(MendelianAlleleProbabilityNHNDeNovoTest.class);
    suite.addTestSuite(SegregationScoreTest.class);
    suite.addTestSuite(SegregationTrivialTest.class);
    suite.addTestSuite(SFunctionTest.class);
    suite.addTestSuite(UniqueIdTest.class);
    return suite;
  }

  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }
}
