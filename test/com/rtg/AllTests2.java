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
package com.rtg;

import junit.framework.Test;
import junit.framework.TestSuite;


/**
 * This class runs the second half of the SLIM tests (the p-z packages).
 *
 */
public class AllTests2 extends TestSuite {

  /**
   * @return suite.
   */
  public static Test suite() {
    final TestSuite suite = new TestSuite("com.rtg");
    suite.addTest(com.rtg.graph.AllTests.suite());
    suite.addTest(com.rtg.position.AllTests.suite());
    suite.addTest(com.rtg.protein.AllTests.suite());
    suite.addTest(com.rtg.reader.AllTests.suite());
    suite.addTest(com.rtg.reference.AllTests.suite());
    suite.addTest(com.rtg.relation.AllTests.suite());
    suite.addTest(com.rtg.report.AllTests.suite());
    suite.addTest(com.rtg.sam.AllTests.suite());
    suite.addTest(com.rtg.visualization.AllTests.suite());
    suite.addTest(com.rtg.scheduler.AllTests.suite());
    suite.addTest(com.rtg.segregation.AllTests.suite());
    suite.addTest(com.rtg.simulation.AllTests.suite());
    suite.addTest(com.rtg.tabix.AllTests.suite());
    suite.addTest(com.rtg.taxonomy.AllTests.suite());
    suite.addTest(com.rtg.util.AllTests.suite());
    suite.addTest(com.rtg.jmx.AllTests.suite());
    suite.addTest(com.rtg.usage.AllTests.suite());
    suite.addTest(com.rtg.variant.AllTests.suite());
    suite.addTest(com.rtg.vcf.AllTests.suite());
    suite.addTest(com.rtg.zooma.AllTests.suite());
    // Uncomment the following lines if you want to see when these tests finish.
    //    suite.addTest(new Test() {
    //
    //      @Override
    //      public int countTestCases() {
    //        return 0;
    //      }
    //
    //      @Override
    //      public void run(TestResult arg0) {
    //        System.out.println("Done2");
    //      }});
    return suite;
  }

  /**
   * Run suite of tests. Used for parallel running.
   * @param args command line arguments - ignored.
   */
  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }
}

