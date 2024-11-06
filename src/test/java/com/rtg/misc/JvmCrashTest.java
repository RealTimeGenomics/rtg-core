/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
package com.rtg.misc;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

import junit.framework.TestCase;

/**
 * Crashed the JVM and shows *** glibc detected *** free(): invalid pointer: 0x000000004022b040 ***
 *
 */
public class JvmCrashTest extends TestCase {

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


