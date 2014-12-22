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
package com.rtg.variant.avr;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 */
public class AllTests extends TestSuite {

  public static Test suite() {
    final TestSuite suite = new TestSuite("com.rtg.variant.avr");
    suite.addTestSuite(AnnotationLoaderTest.class);
    suite.addTestSuite(AttributeExtractorTest.class);
    suite.addTestSuite(BuilderCliTest.class);
    suite.addTestSuite(BooleanConverterTest.class);
    suite.addTestSuite(DoubleConverterTest.class);
    suite.addTestSuite(DerivedAnnotationTest.class);
    suite.addTestSuite(DummyModelBuilderTest.class);
    suite.addTestSuite(DummyPredictModelTest.class);
    suite.addTestSuite(FormatAnnotationTest.class);
    suite.addTestSuite(GtQualComplexMultiplierModelTest.class);
    suite.addTestSuite(GtQualComplexMultiplierModelBuilderTest.class);
    suite.addTestSuite(InfoAnnotationTest.class);
    suite.addTestSuite(IntegerConverterTest.class);
    suite.addTestSuite(MlAvrModelBuilderTest.class);
    suite.addTestSuite(MlAvrPredictModelTest.class);
    suite.addTestSuite(ModelFactoryTest.class);
    suite.addTestSuite(ModelTypeTest.class);
    suite.addTestSuite(NullModelTest.class);
    suite.addTestSuite(NullModelBuilderTest.class);
    suite.addTestSuite(PredictCliTest.class);
    suite.addTestSuite(AvrStatsCliTest.class);
    suite.addTestSuite(QualAnnotationTest.class);
    suite.addTestSuite(BooleanConverterTest.class);
    suite.addTestSuite(DoubleConverterTest.class);
    suite.addTestSuite(IntegerConverterTest.class);
    suite.addTestSuite(VcfDatasetTest.class);
    suite.addTestSuite(StringConverterTest.class);
    return suite;
  }

  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }

}
