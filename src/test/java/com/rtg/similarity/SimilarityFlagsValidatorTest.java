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
package com.rtg.similarity;

import java.io.File;
import java.io.IOException;

import com.rtg.util.TestUtils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.diagnostic.CliDiagnosticListener;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public class SimilarityFlagsValidatorTest extends TestCase {

  public void testValidator() throws IOException {
    Diagnostic.setLogStream();
    final MemoryPrintStream err = new MemoryPrintStream();
    final CliDiagnosticListener listener = new CliDiagnosticListener(err.printStream());
    Diagnostic.addListener(listener);
    final File tempDir = FileHelper.createTempDirectory();
    try {
      final CFlags flags = new CFlags("PhyloTest", TestUtils.getNullPrintStream(), err.printStream());
      SimilarityCli.initFlags(flags);
      checkErrorMessage(flags, new String[] {"-o", "blah", "-I", "ba", "-i", "humbug"}, err, "Only set one of --input or --input-list-file");
      checkErrorMessage(flags, new String[] {"-o", "blah", "-i", "humbug"}, err, "The specified SDF, \"humbug\", does not exist.");
      final File fakePaired = new File(tempDir, "fakePaired");
      assertTrue(fakePaired.mkdir());
      final File left = new File(fakePaired, "left");
      assertTrue(left.mkdir());
      assertTrue(new File(fakePaired, "right").mkdir());
      checkErrorMessage(flags, new String[] {"-o", "blah", "-i", fakePaired.getPath()}, err, "The specified SDF, \"" + fakePaired.getPath() + "\", is a paired end SDF.");
      checkErrorMessage(flags, new String[] {"-o", "blah"}, err, "Must set one of --input or --input-list-file");
      checkErrorMessage(flags, new String[] {"-o", "blah", "-I", "ba"}, err, "The specified list file, \"ba\", does not exist.");
      checkErrorMessage(flags, new String[] {"-o", "blah", "-I", fakePaired.getPath()}, err, "The specified list file,", "\"" + fakePaired.getPath() + "\",", "directory.");
      final File fakeList = new File(tempDir, "fakeList.txt");
      assertTrue(fakeList.createNewFile());
      checkErrorMessage(flags, new String[] {"-o", fakePaired.getPath(), "-I", fakeList.getPath()}, err, "The directory", "\"" + fakePaired.getPath() + "\"", "already exists.");
      checkErrorMessage(flags, new String[] {"-o", "blah", "-I", fakeList.getPath(), "-w", "-1"}, err, "--word must be in the range [1");
      checkErrorMessage(flags, new String[] {"-o", "blah", "-I", fakeList.getPath(), "-w", "0"}, err, "--word must be in the range [1");
      checkErrorMessage(flags, new String[] {"-o", "blah", "-I", fakeList.getPath(), "-w", "33"}, err, "--word must be in the range [1, 32]");
      checkErrorMessage(flags, new String[] {"-o", "blah", "-I", fakeList.getPath(), "--max-reads", "0"}, err, "--max-reads must be greater than 0");
      checkErrorMessage(flags, new String[] {"-o", "blah", "-I", fakeList.getPath(), "-w", "20", "-s", "20"}, err);
      checkErrorMessage(flags, new String[] {"-o", "blah", "-i", left.getPath(), "--max-reads", "1"}, err, "Only set --max-reads when using --input-list-file");
      checkErrorMessage(flags, new String[] {"-o", "blah", "-i", left.getPath()}, err);
    } finally {
      Diagnostic.removeListener(listener);
      err.close();
      assertTrue(FileHelper.deleteAll(tempDir));
    }
  }

  private void checkErrorMessage(CFlags flags, String[] args, MemoryPrintStream err, String... errMessages) {
    err.outputStream().reset();
    if (errMessages.length == 0) {
      assertTrue(flags.setFlags(args));
      assertEquals("", err.toString());
    } else {
      assertFalse(flags.setFlags(args));
      TestUtils.containsAll(err.toString(), errMessages);
    }
  }

}
