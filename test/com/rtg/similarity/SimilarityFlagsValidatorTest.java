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
    CliDiagnosticListener listener = new CliDiagnosticListener(err.printStream());
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
      checkErrorMessage(flags, new String[] {"-o", "blah", "-I", fakeList.getPath(), "-w", "-1"}, err, "The specified flag \"--word\" has invalid value \"-1\". It should be greater than or equal to \"1\".");
      checkErrorMessage(flags, new String[] {"-o", "blah", "-I", fakeList.getPath(), "-w", "0"}, err, "The specified flag \"--word\" has invalid value \"0\". It should be greater than or equal to \"1\".");
      checkErrorMessage(flags, new String[] {"-o", "blah", "-I", fakeList.getPath(), "-w", "33"}, err, "The specified flag \"--word\" has invalid value \"33\". It should be less than or equal to \"32\".");
      checkErrorMessage(flags, new String[] {"-o", "blah", "-I", fakeList.getPath(), "--max-reads", "0"}, err, "The --max-reads must be greater than 0");
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
