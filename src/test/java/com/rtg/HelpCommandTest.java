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
package com.rtg;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Locale;

import com.rtg.launcher.globals.GlobalFlags;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;

import junit.framework.TestCase;

/**
 * Test for corresponding class
 */
public class HelpCommandTest extends TestCase {


  private void basicHelpTester(String input) throws Exception {
    final ByteArrayOutputStream baos = new ByteArrayOutputStream();
    final ByteArrayOutputStream berr = new ByteArrayOutputStream();
    final PrintStream err = new PrintStream(berr);

    try {
      assertEquals(0, HelpCommand.mainInit(new String[]{input}, baos, err, CoreCommand.INFO));
    } finally {
      baos.close();
    }
    berr.flush();
    err.flush();
    final String output = baos.toString();

    for (Command mod : CoreCommand.INFO.commands()) {
      if (mod.isHidden()) {
        assertFalse("Did contain hidden module " + mod.toString().toLowerCase(Locale.getDefault()), output.contains("\t" + mod.toString().toLowerCase(Locale.getDefault()) + " "));
        continue;
      }
      assertTrue("Did not contain " + mod.getCategory().toString(), output.contains(mod.getCategory().getLabel()));
      if (CoreCommand.HELP != mod) {
        assertTrue("Did not contain " + mod.toString().toLowerCase(Locale.getDefault()), output.contains("\t" + mod.toString().toLowerCase(Locale.getDefault())));
      }
    }
  }

  public void testHelp() throws Exception {
    basicHelpTester("help");
  }

  public void testInvalid() throws Exception {
    basicHelpTester("foo");
  }

  public void testOutputUnlicensed() {
    final PrintStream originalErr = System.err;
    final ByteArrayOutputStream bos = new ByteArrayOutputStream();
    final PrintStream ps = new PrintStream(bos);
    System.setErr(ps);
    try {
      HelpCommand.outputUnlicensedModule(ToolsCommand.GENOMESIM);
      fail();
    } catch (RuntimeException e) {
      //expected
    }
    ps.flush();
    final String output = bos.toString();
    ps.close();
    final String expectedOutput = "The GENOMESIM command has not been enabled by your current license." + StringUtils.LS
                                + "Please contact support@realtimegenomics.com to have this command licensed." + StringUtils.LS;
    assertEquals(expectedOutput, output);
    System.setErr(originalErr);
  }

  public void testNull() throws Exception {
    final ByteArrayOutputStream baos = new ByteArrayOutputStream();
    final ByteArrayOutputStream berr = new ByteArrayOutputStream();
    final PrintStream err = new PrintStream(berr);
    try {
      GlobalFlags.resetAccessedStatus();
      assertEquals(1, HelpCommand.mainInit(null, baos, err, CoreCommand.INFO));
    } finally {
      baos.close();
    }
    berr.flush();
    err.flush();
    assertEquals(HelpCommand.getUsage(false, CoreCommand.INFO), berr.toString());
    assertEquals("", baos.toString());
  }

  public void testValid() throws Exception {
    final ByteArrayOutputStream baos = new ByteArrayOutputStream();
    final ByteArrayOutputStream berr = new ByteArrayOutputStream();
    final PrintStream err = new PrintStream(berr);
    try {
      GlobalFlags.resetAccessedStatus();
      assertEquals(0, HelpCommand.mainInit(new String[0], baos, err, CoreCommand.INFO));
    } finally {
      baos.close();
    }
    final String output = baos.toString();
    TestUtils.containsAll(output, StringUtils.LS + StringUtils.LS + "Data formatting:"
                                , StringUtils.LS + StringUtils.LS + "Read mapping:"
                                , StringUtils.LS + StringUtils.LS + "Variant detection:"
                                , StringUtils.LS + StringUtils.LS + "Metagenomics:"
                                , StringUtils.LS + StringUtils.LS + "Simulation:"
                                , StringUtils.LS + StringUtils.LS + "Utility:");

    final Command[] arr = CoreCommand.INFO.commands();
    String shortestName = arr[0].toString();
    for (Command mod : arr) {
      if (mod.isHidden()) {
        continue;
      }
      if (mod.toString().length() < shortestName.length()) {
        shortestName = mod.toString();
      }
      assertTrue("Did not contain " + mod.getCategory().toString(), output.contains(mod.getCategory().getLabel()));
      if (CoreCommand.HELP != mod) {
        assertTrue("Did not contain " + mod.toString().toLowerCase(Locale.getDefault()), output.contains(mod.toString().toLowerCase(Locale.getDefault())));
      }
    }
  }

  public void testModuleHelp() throws IOException {
    final ByteArrayOutputStream out = new ByteArrayOutputStream();
    final ByteArrayOutputStream berr = new ByteArrayOutputStream();
    final PrintStream err = new PrintStream(berr);
    GlobalFlags.resetAccessedStatus();
    ToolsCommand.FORMAT.mainInit(new String[]{"-h"}, out, err);

    final ByteArrayOutputStream bhelp = new ByteArrayOutputStream();
    final ByteArrayOutputStream bhelpout = new ByteArrayOutputStream();
    try (PrintStream help = new PrintStream(bhelp)) {
      GlobalFlags.resetAccessedStatus();
      assertEquals(0, HelpCommand.mainInit(new String[]{"format"}, bhelpout, help, CoreCommand.INFO));
      help.flush();
    } finally {
      bhelp.close();
      bhelpout.close();
    }
    assertEquals(bhelpout.toString(), out.toString());
  }

  private void xHelpTester(String xflag) throws IOException {
    final ByteArrayOutputStream baos = new ByteArrayOutputStream();
    final ByteArrayOutputStream berr = new ByteArrayOutputStream();
    final PrintStream err = new PrintStream(berr);
    try {
      GlobalFlags.resetAccessedStatus();
      assertEquals(0, HelpCommand.mainInit(new String[]{xflag}, baos, err, CoreCommand.INFO));
    } finally {
      baos.close();
    }
    final String output = baos.toString();

    TestUtils.containsAll(output, "Usage: rtg COMMAND [OPTION]..."
                                , "Type 'rtg help COMMAND' for help on a specific command."
                                , "The following commands are available:");

    String longestName = "";
    String desc = "";
    for (Command mod : CoreCommand.INFO.commands()) {
      if (mod.toString().length() > longestName.length()) {
        longestName = mod.toString();
        desc = mod.getCommandDescription();
      }
      assertTrue("Did not contain " + mod.getCategory().toString(), output.contains(mod.getCategory().getLabel()));
      if (CoreCommand.HELP != mod) {
        assertTrue("Did not contain " + mod.toString().toLowerCase(Locale.getDefault()), output.contains(mod.toString().toLowerCase(Locale.getDefault())));
      }
    }
    final String exp = longestName.toLowerCase(Locale.getDefault()) + " " + desc;
    final String out = output.replaceAll(StringUtils.LS,  "").replaceAll("\\s+", " ");
    assertTrue(exp + "\nNot in \n" + out, out.contains(exp));
  }

  public void testValidHidden() throws Exception {
    xHelpTester("--Xhidden");
  }

  public void testValidXHelp() throws Exception {
    xHelpTester("--Xhelp");
  }

  public void testGetUsageString() {
    final String output = HelpCommand.getUsage(false, 60, CoreCommand.INFO);
    TestUtils.containsAll(output, "\tmap ", "\tread mapping");
  }

  public void testLongestLengthModule() {
    assertEquals(0, HelpCommand.getLongestLengthModule(new Command[0]));
    assertEquals(6, HelpCommand.getLongestLengthModule(new Command[]{ToolsCommand.FORMAT}));
  }
}
