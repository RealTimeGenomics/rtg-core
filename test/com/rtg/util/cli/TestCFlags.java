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
package com.rtg.util.cli;

import java.io.IOException;
import java.util.HashSet;
import java.util.Locale;

import com.rtg.util.TestUtils;

import junit.framework.Assert;

/**
 * Test a <code>CFlags</code> object for various metrics used in SLIM.
 *
 */
public final class TestCFlags {

  private TestCFlags() { }

  private static final String[] OUTLAWS = {
    "preread",    // should be SDF
    "pre-read",   // should be SDF
    "specifies",  // probably redundant
    "don't",      // colloquial
    "CG",         // should be Complete Genomics
    "classname",  // users don't know about Java
    "class",      // users don't know about Java
    "licence",    // we are an American company
    "colour",     // we are an American company
  };

  private static final HashSet<String> PDESC = new HashSet<>();
  static {
    PDESC.add(null);
    PDESC.add("BOOL");
    PDESC.add("DIR");
    PDESC.add("EXPRESSION");
    PDESC.add("FILE");
    PDESC.add("FLOAT");
    PDESC.add("FORMAT");
    PDESC.add("INT");
    PDESC.add("MODEL");
    PDESC.add("NAME");
    PDESC.add("STRING");
    PDESC.add("PERCENTAGE");
    PDESC.add("SDF");
    PDESC.add("SDF|FILE");
    PDESC.add("SEX");
    PDESC.add("STRING|FILE");
  }

  private static void checkDescriptionConstraints(final Flag f) {
    // Most of these conventions are based on what is seen in standard Unix commands
    // and man pages
    final String desc = f.getDescription();
    //System.err.println("Description: " + desc);
    final String name = f.getName();
    Assert.assertNotNull("Null description: --" + name, desc);
    if ((name != null) && (name.charAt(0) == 'X')) {
      return;
    }
    Assert.assertTrue("Description is too short: --" + name + " desc: " + desc, desc.length() > 8);
    Assert.assertTrue("Description should start with lowercase: --" + name + " desc: " + desc, Character.isLowerCase(desc.charAt(0)) || Character.isUpperCase(desc.charAt(1)));
    Assert.assertTrue("Description should start with lowercase: --" + name + " desc: " + desc, Character.isLetterOrDigit(desc.charAt(0)));
    Assert.assertFalse("Description should not end with \".\": --" + name + " desc: " + desc, desc.endsWith("."));
    final String[] parts = desc.split("\\s.,;:\"!?");
    for (final String o : OUTLAWS) {
      for (final String p : parts) {
        Assert.assertFalse("Outlawed \"" + o + "\" occurs in flag --" + name + " usage: " + desc, o.equals(p));
      }
    }
    CheckSpelling.check(name, desc);
    final String lcDesc = desc.toLowerCase(Locale.getDefault());

    // Should not mention the word default twice
    final int d = lcDesc.indexOf("default");
    if (d != -1) {
      Assert.assertEquals(d, lcDesc.lastIndexOf("default"));
    }
    if (name != null) {
      if (name.length() < 3) {
        Assert.fail("Long flag name is too short: --" + name);
      }
      for (int k = 0; k < name.length(); k++) {
        if (Character.isWhitespace(name.charAt(k))) {
          Assert.fail("Name of flag contains whitespace: --" + name);
        }
      }
    }
    final String pd = f.getParameterDescription();
    Assert.assertTrue(pd + " is not an allowed parameter description", PDESC.contains(pd));
  }

  /**
   * Check various syntactic properties of a <code>CFlags</code> description.
   *
   * @param flags the flags
   * @param contains strings required to be present
   */
  public static void check(final CFlags flags, final String... contains) {
    Assert.assertNotNull(flags);
    final StringBuilder problems = new StringBuilder();
    try {
      CheckSpelling.setSpelling(problems);
      for (final Flag f : flags.getRequired()) {
        checkDescriptionConstraints(f);
      }
      for (final Flag f : flags.getOptional()) {
        checkDescriptionConstraints(f);
      }
      if (problems.length() > 0) {
        Assert.fail(problems.toString());
      }
    } catch (final IOException e) {
      Assert.fail(e.getMessage());
    }
    //System.err.println("Usage TestCFlags.java\n" + flags.getUsageString());
    final String usage = flags.getUsageString().replaceAll("    --", "\\\\0, --").replaceAll("\\s+", " ");
    Assert.assertNotNull(usage);
    if (contains != null) {
      TestUtils.containsAll(usage, contains);
    }
  }

  /**
   * Check various syntactic properties of a <code>CFlags</code> description.
   *
   * @param flags the flags
   * @param contains strings required to be present
   */
  public static void checkExtendedUsage(final CFlags flags, final String... contains) {
    Assert.assertNotNull(flags);
    //System.err.println("Usage TestCFlags.java\n" + flags.getUsageString());
    final String usage = flags.getExtendedUsageString().replaceAll("\\s+", " ");
    Assert.assertNotNull(usage);
    if (contains != null) {
      TestUtils.containsAll(usage, contains);
    }
  }

}

