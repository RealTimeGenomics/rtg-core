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
package com.rtg.util.io;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 * Test class for all tests in this directory.
 *
 */
public class AllTests extends TestSuite {

  public static Test suite() {
    final TestSuite suite = new TestSuite("com.rtg.util.io");
    suite.addTestSuite(ConcurrentByteQueueTest.class);
    suite.addTestSuite(BufferedRandomAccessFileTest.class);
    suite.addTestSuite(BufferedOutputStreamFixTest.class);
    suite.addTestSuite(AsynchInputStreamTest.class);
    suite.addTestSuite(AsynchOutputStreamTest.class);
    suite.addTestSuite(ClosedFileInputStreamTest.class);
    suite.addTestSuite(GzipAsynchInputStreamTest.class);
    suite.addTestSuite(GzipAsynchOutputStreamTest.class);
    suite.addTestSuite(IOUtilsTest.class);
    suite.addTestSuite(LogFileTest.class);
    suite.addTestSuite(LineWriterTest.class);
    suite.addTestSuite(LogSimpleTest.class);
    suite.addTestSuite(ByteArrayIOUtilsTest.class);
    suite.addTestSuite(FileUtilsTest.class);
    suite.addTestSuite(MemoryPrintStreamTest.class);
    suite.addTestSuite(InputFileUtilsTest.class);
    suite.addTestSuite(PartitionTest.class);
    suite.addTestSuite(SimpleArchiveTest.class);
    suite.addTestSuite(AdjustableGZIPOutputStreamTest.class);
    suite.addTestSuite(RandomAccessFileStreamTest.class);
    suite.addTestSuite(FalseSeekableStreamTest.class);
    suite.addTestSuite(TestDirectoryTest.class);
    suite.addTest(com.rtg.util.io.bzip2.AllTests.suite());
    return suite;
  }
}

