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

package com.rtg.visualization;

import java.io.File;
import java.io.IOException;

import com.rtg.util.cli.CFlags;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;


/**
 */
public class AviewParamsTest extends TestCase {

  public void testAviewParams() throws IOException {
    final File f = FileUtils.createTempDir("aview", "params");
    try {
      final File foo = File.createTempFile("foo", "def", f);
      final File bar = File.createTempFile("bar", "def", f);
    final String[] args = {
      foo.getPath(), bar.getPath(), "-r", "read1", "-r", "read2",
      "-c", "snps.gz", "-b", "generated.snps", "--region", "really long name:1000+1000", "-t", "template",
      "--print-reference-line", "100", "-U", "unmapped1", "-U", "unmapped2"
    };
    final CFlags flags = new CFlags(); //("AviewParamsTest", System.out, System.err);
    AviewParams.initFlags(flags);
    assertTrue(flags.setFlags(args));
    final AviewParams parameters = AviewParams.makeParams(flags);
    assertEquals(2, parameters.alignmentsFiles().length);
    assertEquals(foo.getPath(), parameters.alignmentsFiles()[0].toString());
    assertEquals(bar.getPath(), parameters.alignmentsFiles()[1].toString());

    assertEquals(2, parameters.unmappedFiles().length);
    assertEquals("unmapped1", parameters.unmappedFiles()[0].getName());
    assertEquals("unmapped2", parameters.unmappedFiles()[1].getName());

    assertEquals(2, parameters.readFiles().length);
    assertEquals("read1", parameters.readFiles()[0].getName());
    assertEquals("read2", parameters.readFiles()[1].getName());

    assertEquals(1000, parameters.start());
    assertEquals(2000, parameters.end());

    assertEquals("snps.gz", parameters.trackFiles()[0].getName());
    assertEquals("generated.snps", parameters.baselineFile().getName());

    assertEquals("really long name", parameters.sequenceName());
    assertEquals("template", parameters.referenceFile().getName());
    assertTrue(parameters.displayDots());
    assertEquals(100, parameters.headerLineRepeat());
    } finally {
      assertTrue(FileHelper.deleteAll(f));
    }
  }

  public void testFullAviewParams() throws IOException {
    final File f = FileUtils.createTempDir("aview", "params");
    try {
      final File foo = File.createTempFile("foo", "def", f);
      final File bar = File.createTempFile("bar", "def", f);
    final String[] args = {
      foo.getPath(), bar.getPath(), "--reads", "read1", "--reads", "read2",
      "--calls", "snps.gz", "--baseline", "generated.snps",
      "--region", "really long name:1000+1000", "--template", "template", "--print-reference-line", "200",
      "--no-color", "--Xmax-as-mated", "10", "--Xmax-as-unmated", "8", "--Xmax-hits", "2", "--no-dots", "--print-cigars",
      "--print-names",
      "--Xunmapped", "unmapped1", "--Xunmapped", "unmapped2"
    };
    final CFlags flags = new CFlags(); //("AviewParamsTest", System.out, System.err);
    AviewParams.initFlags(flags);
    assertTrue(flags.setFlags(args));
    final AviewParams p = AviewParams.makeParams(flags);
    assertEquals(2, p.alignmentsFiles().length);
    assertEquals(foo.getPath(), p.alignmentsFiles()[0].toString());
    assertEquals(bar.getPath(), p.alignmentsFiles()[1].toString());

    assertEquals(2, p.readFiles().length);
    assertEquals("read1", p.readFiles()[0].getName());
    assertEquals("read2", p.readFiles()[1].getName());

    assertEquals(2, p.unmappedFiles().length);
    assertEquals("unmapped1", p.unmappedFiles()[0].getName());
    assertEquals("unmapped2", p.unmappedFiles()[1].getName());

    assertEquals(1000, p.start());
    assertEquals(2000, p.end());

    assertEquals("snps.gz", p.trackFiles()[0].getName());
    assertEquals("generated.snps", p.baselineFile().getName());

    assertEquals("really long name", p.sequenceName());
    assertEquals("template", p.referenceFile().getName());

    assertEquals(200, p.headerLineRepeat());
    assertFalse(p.displayDots());

    assertFalse(p.useTerminalColor());
    assertTrue(p.printCigars());
    assertEquals(10, p.maxMatedAlignmentScore());
    assertEquals(8, p.maxUnmatedAlignmentScore());
    assertEquals(2, p.maxIhScore());
    assertTrue(p.printReadName());
    } finally {
      assertTrue(FileHelper.deleteAll(f));
    }
  }

  public void testDefaultParams() throws IOException {
    final File f = FileUtils.createTempDir("aview", "params");
    try {
      final File foo = File.createTempFile("foo", "def", f);
    final String[] arguments = {
        foo.getPath(), "-c", "snps.gz", "--region", "really long name:1000+1000", "-t", "template",
      };
      final CFlags flags = new CFlags(); //("AviewParamsTest", System.out, System.err);
      AviewParams.initFlags(flags);
      assertTrue(flags.setFlags(arguments));
      final AviewParams params = AviewParams.makeParams(flags);
      assertEquals(1, params.alignmentsFiles().length);
      assertEquals(foo.getPath(), params.alignmentsFiles()[0].getPath());

      assertEquals(0, params.readFiles().length);
      assertEquals(0, params.unmappedFiles().length);

      assertEquals(1000, params.start());
      assertEquals(2000, params.end());

      assertEquals("snps.gz", params.trackFiles()[0].getName());
      assertNull(params.baselineFile());
      assertEquals("really long name", params.sequenceName());
      assertEquals("template", params.referenceFile().getName());
      assertTrue(params.displayDots());
      assertEquals(0, params.headerLineRepeat());
      assertTrue(params.useTerminalColor());
      assertFalse(params.printCigars());
      assertEquals(Integer.MAX_VALUE, params.maxMatedAlignmentScore());
      assertEquals(Integer.MAX_VALUE, params.maxUnmatedAlignmentScore());
      assertEquals(Integer.MAX_VALUE, params.maxIhScore());
      assertFalse(params.printReadName());

      assertNotNull(new AviewParamsBuilder().create());
    } finally {
      assertTrue(FileHelper.deleteAll(f));
    }
  }

}
