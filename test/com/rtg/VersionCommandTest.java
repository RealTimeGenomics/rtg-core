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

import com.rtg.util.Constants;
import com.rtg.util.Environment;
import com.rtg.util.License;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.io.MemoryPrintStream;

import junit.framework.TestCase;

/**
 * Tests for the version command.
 */
public class VersionCommandTest extends TestCase {

  public void ramStringTest(final String string) {
    final Double actualpercentage = Runtime.getRuntime().maxMemory() * 100 / (Environment.getTotalMemory() * 1.0d);
    final double maxMemory = Runtime.getRuntime().maxMemory() * 10.0 / 1024 / 1024 / 1024;
    assertTrue(string.contains("RAM: "));
    assertTrue(string.contains("GB RAM can be used by rtg"));
    assertTrue(string.contains("RAM: " + ((int) maxMemory / 10.0) + "GB of"));

    final String[] lines = string.split("\n");
    for (final String line : lines) {
      if (line.contains("RAM: ")) {
        final String[] chunks = line.trim().split(" ");
        assertTrue(chunks[3].endsWith("GB"));
        final int maxpos = chunks[3].indexOf("GB");
        final String max = chunks[3].substring(0, maxpos);
        final Double maxmem = Double.parseDouble(max);
        assertTrue(maxmem >= 0);
        assertTrue(chunks[1].endsWith("GB"));
        final int gbpos = chunks[1].indexOf("GB");
        final String alloc = chunks[1].substring(0, gbpos);
        final Double allocmem = Double.parseDouble(alloc);
        assertTrue(allocmem >= 0);
        assertTrue(allocmem <= maxmem);


        assertTrue(chunks[10].startsWith("("));
        assertTrue(chunks[10].endsWith("%)"));
        final String percent = chunks[10].substring(1, chunks[10].length() - 2);
        final Double percentdbl = Double.parseDouble(percent);
        assertEquals(actualpercentage.intValue(), percentdbl.intValue());
      }
    }
    String ramString = VersionCommand.getRamString(100, 1000, 1000);
    assertEquals("100.0GB of 100.0GB RAM can be used by rtg (100%)", ramString);
    ramString = VersionCommand.getRamString(9, 90, 1000);
    assertEquals("9.0GB of 100.0GB RAM can be used by rtg (9%)", ramString);
    ramString = VersionCommand.getRamString(80, 80, 100);
    assertEquals("8.0GB of 10.0GB RAM can be used by rtg (80%)", ramString);
  }

  public void testStrings() {
    final MemoryPrintStream bout = new MemoryPrintStream();
    try {
      assertEquals(0, VersionCommand.mainInit(bout.printStream()));
      bout.outputStream().reset();
      assertEquals(0, VersionCommand.mainInit(new String[]{}, bout.outputStream()));
    } finally {
      bout.close();
    }

    final String versionOut = bout.toString();
    assertTrue(versionOut.startsWith("Product: " + Environment.getProductName() + StringUtils.LS
        + "Core Version: " + Environment.getCoreVersion() + StringUtils.LS));
    ramStringTest(versionOut);

    TestUtils.containsAll(versionOut, "License: " + License.getMessage() + StringUtils.LS,
        "JVM: ",
        "Contact: " + Constants.SUPPORT_EMAIL_ADDR + StringUtils.LS,
        "Patents / Patents pending:",
        "Citation:",
        "(c) Real Time Genomics, 201");

  }

}
