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

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.net.InetAddress;
import java.net.UnknownHostException;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;
import java.util.Map;
import java.util.ResourceBundle;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;


public class EnvironmentTest extends TestCase {

  public EnvironmentTest(final String name) {
    super(name);
  }

  public static Test suite() {
    return new TestSuite(EnvironmentTest.class);
  }
  /**
   * Main to run from tests from command line.
   * @param args ignored.
   */
  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }

  public void test() {
    HashSet<String> emptyKeys = new HashSet<>();
    emptyKeys.add("user.variant");
    emptyKeys.add("exe4j.tempDir");
    emptyKeys.add("user.timezone");
    emptyKeys.add("slim.bl2seq");
    emptyKeys.add("sun.cpu.isalist");
    emptyKeys.add("sun.os.patch.level");
    emptyKeys.add("user.script");
    emptyKeys.add("regression.update");

    String vmName = System.getProperty("java.vm.name", "");
    if (vmName.startsWith("IBM J9")) {
      emptyKeys.add("com.ibm.util.extralibs.properties");
      emptyKeys.add("java.awt.fonts");
      emptyKeys.add("com.ibm.jcl.checkclasspath");
      emptyKeys.add("sun.java2d.fontpath");
    }

    HashMap<String, String> hs = Environment.getEnvironmentMap();
    //Set<String> st = hs.keySet();
    //String[] keys = st.toArray(new String[st.size()]);
    for (final Map.Entry<String, String> entry : hs.entrySet()) {
      final String s = entry.getKey();
      //System.err.println(s + ":" + hs.get(s));
      if (!emptyKeys.contains(s.toLowerCase(Locale.getDefault()))) {
        assertTrue(s, entry.getValue().length() > 0);
      }
    }
  }

  public void testkeys() {

    assertEquals(Environment.RUNTIME_FREEMEM, "runtime.freememory");
    assertEquals(Environment.RUNTIME_MAXMEM, "runtime.maxmemory");
    assertEquals(Environment.RUNTIME_TOTALMEM, "runtime.totalMemory");

    assertEquals(Environment.PROCESSOR_ARCH, "processor.arch");
    assertEquals(Environment.PROCESSOR_COUNT, "processor.count");

    assertEquals(Environment.OS_FREEMEM, "os.freememory");
    assertEquals(Environment.OS_FREEMEM, "os.freememory");
    assertEquals(Environment.OS_TOTALMEM, "os.totalmemory");

    HashMap<String, String> hs = Environment.getEnvironmentMap();
    Collection<String> st = hs.keySet();

    assertTrue(!st.isEmpty());
    assertTrue(!hs.isEmpty());
    assertTrue(hs.size() > 8);
    assertTrue(hs.containsKey("runtime.freememory"));
    assertTrue(hs.containsKey("runtime.maxmemory"));

    assertTrue(hs.containsKey("runtime.totalMemory"));
    assertTrue(hs.containsKey("processor.arch"));
    /*
    //Test commented out due to its inabilty to compile on anything other than the Sun JDK
    final OperatingSystemMXBean bean = ManagementFactory.getOperatingSystemMXBean();
    if (bean instanceof com.sun.management.OperatingSystemMXBean) {
      assertTrue(hs.containsKey(Environment.OS_FREEMEM));
      assertTrue(hs.containsKey("processor.count"));
      // SAI comment out buggy test, memory can change in between the calls!
      //      assertEquals(Environment.getFreeMemory(), ((com.sun.management.OperatingSystemMXBean) bean).getFreePhysicalMemorySize());

      assertEquals(Environment.getTotalMemory(), ((com.sun.management.OperatingSystemMXBean) bean).getTotalPhysicalMemorySize());
      assertEquals(Environment.getAvailableProcessors(), bean.getAvailableProcessors());
     }
    */

    assertTrue(Environment.getTotalMemory() > 500000000); //more than 500M
  }

  public void testVaules() {
    HashMap<String, String> hs = Environment.getEnvironmentMap();
    hs.put(Environment.OS_FREEMEM, "200000");
    assertEquals(hs.get(Environment.OS_FREEMEM), "200000");
  }

  public void testMain0() throws Exception {
    final PrintStream oldOut = System.out;
    try {
      final ByteArrayOutputStream bos = new ByteArrayOutputStream();
      try {
        try (PrintStream out = new PrintStream(bos)) {
          System.setOut(out);
          Environment.main();
        }
      } finally {
        bos.close();
      }
      final String s = bos.toString();
      assertTrue(s.contains("file.separator"));
      assertTrue(s.contains("vendor"));
      assertTrue(s.contains("lang.version"));
      assertTrue(s.contains("os.arch"));
      assertTrue(s.contains("os.freememory"));
      assertTrue(s.contains("os.name"));
      assertTrue(s.contains("runtime.freememory"));
      assertTrue(s.contains("runtime.maxmemory"));
      assertTrue(s.contains("runtime.totalMemory"));
      assertTrue(s.contains("user.home"));
    } finally {
      System.setOut(oldOut);
    }
  }


  public void testVersion() {
    final String product = Environment.getVersion("com/rtg/product.name");
    final String svn = Environment.getVersion("com/rtg/slim.version");
    final ResourceBundle environment = ResourceBundle.getBundle("com.rtg.util.Environment", Locale.getDefault());
    String version = environment.getString("VERSION");
    final String time = Environment.getVersion("com/rtg/build.time");
    //System.err.println(Environment.getVersion());
    assertEquals(version, Environment.getRelease());
    assertEquals(product + " / Core " + version + " build " + svn + " (" + time + ")", Environment.getVersion());
    assertFalse(Environment.getVersion().contains("cannot"));

  }


  private static final String HOST_NAME = "host.name";

  public void testHostname() {
    String expected;
    try {
      expected = InetAddress.getLocalHost().getHostName();
    } catch (UnknownHostException ex) {
      expected = "<Unknown Host>";
    }
    assertTrue(Environment.getEnvironmentMap().containsKey(HOST_NAME));
    final String actual = Environment.getEnvironmentMap().get(HOST_NAME);
    assertEquals(expected, actual);
  }
}

