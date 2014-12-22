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
package com.rtg.util;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Utils for spawning a separate JVM.
 *
 */
public final class SpawnJvm {

  // private constructor so no instances can be made
  private SpawnJvm() { }

  /**
   * Spawn a process.  Maximum memory is set to 64 MiB.
   *
   * @param className class to run
   * @param args arguments to <code>main</code>
   * @return the Process object
   * @throws IOException if there is an I/O problem.
   */
  public static Process spawn(final String className, final String... args) throws IOException {
    return spawn(64 * 1025 * 1024L, className, args);
  }

  /**
   * Spawn a process with specified memory allocation.
   *
   * @param maxMem maximum amount of memory to allocate
   * @param className class to run
   * @param args arguments to <code>main</code>
   * @return the Process object
   * @throws IOException if there is an I/O problem.
   */
  public static Process spawn(final long maxMem, final String className, final String... args) throws IOException {
    return spawn(maxMem, className, false, args);
  }

  /**
   * Spawn a process with specified memory allocation.
   *
   * @param maxMem maximum amount of memory to allocate
   * @param className class to run
   * @param enableAssert assertions
   * @param args arguments to <code>main</code>
   * @return the Process object
   * @throws IOException if there is an I/O problem.
   */
  public static Process spawn(final long maxMem, final String className, final boolean enableAssert, final String... args) throws IOException {
    final String slash = System.getProperty("file.separator");
    final String javahome = System.getProperty("java.home");
    final String classpath = System.getProperty("java.class.path");

    final ArrayList<String> command = new ArrayList<>();
    command.add(javahome + slash + "bin" + slash + "java");
    if (enableAssert) {
      command.add("-ea");
    }
//    command.add("-server");
    command.add("-Xmx" + maxMem);
    //command.add("-Xms" + maxMem); //doesn't work on Linux ok on MacOSX back out till we understand what is going on - JC
    command.add("-cp");
    command.add(classpath);
    command.add(className);
    command.addAll(Arrays.asList(args));
    return new ProcessBuilder(command).start();
  }

}


