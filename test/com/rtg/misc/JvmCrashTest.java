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
package com.rtg.misc;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Crashed the JVM and shows *** glibc detected *** free(): invalid pointer: 0x000000004022b040 ***
 *
 */
public class JvmCrashTest extends TestCase {

  public JvmCrashTest(final String name) {
    super(name);
  }

  public static Test suite() {
    return new TestSuite(JvmCrashTest.class);
  }

  /**
   * Main to run from tests from command line.
   * @param args ignored.
   */
  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }

  public void testJvm() throws IOException {
    runClassBothStreams(134217728, "com.rtg.util.SpawnJvmTest$MaxMem");
  }

  public static void runClassBothStreams(final int maxMem, final String className, final String... args) throws IOException {
    spawn(maxMem * 1024L * 1024L, className, args);
  }

  public static Process spawn(final long maxMem, final String className, final String... args) throws IOException {
    final String slash = System.getProperty("file.separator");
    final String javahome = System.getProperty("java.home");
    final String classpath = System.getProperty("java.class.path");

    final ArrayList<String> command = new ArrayList<>();
    command.add(javahome + slash + "bin" + slash + "java");
    command.add("-Xmx" + maxMem);
    // command.add("-Xms" + maxMem); //doesn't work on Linux ok on MacOSX back out till we understand what is going on - JC
    command.add("-cp");
    command.add(classpath);
    command.add(className);
    command.addAll(Arrays.asList(args));
    final ProcessBuilder pb = new ProcessBuilder(command);
    final Process pr = pb.start();
    pr.getErrorStream().close();
    pr.getInputStream().close();
    pr.getOutputStream().close();
    return pr;
  }
}


