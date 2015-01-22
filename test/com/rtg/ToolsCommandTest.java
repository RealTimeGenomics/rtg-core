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

import java.util.Collections;
import java.util.HashSet;

import com.rtg.util.TestUtils;

import junit.framework.TestCase;

/**
 */
public class ToolsCommandTest extends TestCase {
  public void testEnum() {
    TestUtils.testEnum(ToolsCommand.class, "[FORMAT, SDF2FASTA, SDF2FASTQ, BGZIP, INDEX, EXTRACT, SDFSTATS, "
            + "SDFSUBSET, SDFSUBSEQ, MENDELIAN, VCFSTATS, VCFMERGE, VCFSUBSET, VCFFILTER, VCFANNOTATE, "
            + "VCFEVAL, PEDFILTER, PEDSTATS, ROCPLOT, VERSION, LICENSE, HELP]");
    final HashSet<Command> displayCommands = new HashSet<>();
    Collections.addAll(displayCommands, ToolsCommand.INFO.commands());
    for (ToolsCommand mod : ToolsCommand.values()) {
      assertTrue(mod.toString(), displayCommands.contains(mod.module()));
    }
  }
}
