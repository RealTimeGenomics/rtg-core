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
package com.rtg.launcher;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;

import com.rtg.mode.DNAFastaSymbolTable;
import com.rtg.mode.ProteinFastaSymbolTable;
import com.rtg.reader.FastaSequenceDataSource;
import com.rtg.reader.PrereadType;
import com.rtg.reader.SequencesWriter;
import com.rtg.util.Utils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.SlimException;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;
import com.rtg.util.test.MockEventListener;

import junit.framework.TestCase;

/**
 */
public class DefaultReaderParamsTest extends TestCase {

  public static DefaultReaderParams createDefaultReaderParams(final File sequenceDir, LongRange readerRestriction) {
    return new DefaultReaderParams(sequenceDir, readerRestriction, false, false, false);
  }

  public DefaultReaderParamsTest(final String name) {
    super(name);
  }

  private void getReaderDNA(final String inputDnaSequence, final File dir) throws IOException {
    final InputStream fqis = new ByteArrayInputStream(inputDnaSequence.getBytes());
    final FastaSequenceDataSource ds = new FastaSequenceDataSource(fqis, new DNAFastaSymbolTable());
    final SequencesWriter sequenceWriter = new SequencesWriter(ds, dir, 20, PrereadType.UNKNOWN, false);
    sequenceWriter.processSequences();

  }

  private void getReaderProtein(final String inputProteinSequence, final File dir) throws IOException {
    final InputStream fqis = new ByteArrayInputStream(inputProteinSequence.getBytes());
    final FastaSequenceDataSource ds = new FastaSequenceDataSource(fqis, new ProteinFastaSymbolTable());
    final SequencesWriter sw = new SequencesWriter(ds, dir, 20, PrereadType.UNKNOWN, false);
    sw.processSequences();

  }

  public void testConstructor() throws IOException {
    Diagnostic.setLogStream();
    final File f = FileUtils.createTempDir("test", "defreadparampro");
    try {
      final File f2 = FileUtils.createTempDir("twst", "defreadparamdna");
      try {
        getReaderProtein(">test\nACNGT\n", f);
        getReaderDNA(">test\nACGT\n", f2);
        final DefaultReaderParams drp = createDefaultReaderParams(f, LongRange.NONE);
        final DefaultReaderParams drp2 = createDefaultReaderParams(f, LongRange.NONE);
        final DefaultReaderParams drp3 = createDefaultReaderParams(f2, LongRange.NONE);
        assertEquals(f, drp.directory());
        assertEquals(5, drp.maxLength());
        assertEquals(Utils.hash(new Object[]{f, Boolean.FALSE}), drp.hashCode());
        assertEquals("ReaderParams directory=" + f + " usememory=" + Boolean.FALSE, drp.toString());
        assertTrue(Arrays.toString(drp.lengths()), Arrays.equals(new int[] {5}, drp.lengths()));
        assertTrue(drp.lengths() == drp.lengths()); //make sure cacheing correctly

        assertTrue(drp.equals(drp));
        assertFalse(drp.equals(null));
        assertTrue(drp.equals(drp2));
        assertFalse(drp.equals(drp3));
        assertTrue(drp.integrity());
        drp.close();
        drp2.close();
        drp3.close();
      } finally {
        FileHelper.deleteAll(f2);
      }
    } finally {
      FileHelper.deleteAll(f);
    }
  }

  public void testDirectoryNotExist() throws IOException {
    final File tempDir = FileUtils.createTempDir("DefaultReaderParamsTest", null);
    try {
      final MockEventListener ev3 = new MockEventListener();
      Diagnostic.addListener(ev3);
      Diagnostic.setLogStream();

      try {
        final File t = new File(tempDir, "t");
        try {
          final DefaultReaderParams rp = createDefaultReaderParams(t, LongRange.NONE);
          rp.reader(); //force creation of reader
          rp.integrity();
          fail();
        } catch (final SlimException ex) {
          ex.printErrorNoLog();
          assertTrue(ev3.compareErrorMessage("Error: The specified SDF, \"" + t.toString() + "\", does not exist."));
        }
      } finally {
        Diagnostic.removeListener(ev3);
      }
    } finally {
      FileHelper.deleteAll(tempDir);
    }
  }

  public void testDirectoryInvalid() throws IOException {
    final MockEventListener ev3 = new MockEventListener();
    Diagnostic.addListener(ev3);
    Diagnostic.setLogStream();
    final File dir = FileUtils.createTempDir("defaultreader", "test");
    try {
      try {
        final DefaultReaderParams rp = createDefaultReaderParams(dir, LongRange.NONE);
        rp.reader(); //force creation of reader
        rp.integrity();
        fail();
      } catch (final SlimException ex) {
        ex.printErrorNoLog();
        assertTrue(ev3.compareErrorMessage("Error: The specified SDF, \"" + dir.toString() + "\", does not seem to contain a valid SDF index."));
      }
    } finally {
      Diagnostic.removeListener(ev3);
    }
    FileHelper.deleteAll(dir);
  }

  public void testDirectoryIsFile() throws IOException {
    final MockEventListener ev3 = new MockEventListener();
    Diagnostic.addListener(ev3);
    Diagnostic.setLogStream();
    final File file = File.createTempFile("defaultreader", "test");
    try {
      try {
        final DefaultReaderParams rp = createDefaultReaderParams(file, LongRange.NONE);
        rp.reader(); //force creation of reader
        rp.integrity();
        fail();
      } catch (final SlimException ex) {
        ex.printErrorNoLog();
        assertTrue(ev3.compareErrorMessage("Error: The specified file, \"" + file.toString() + "\", is not an SDF."));
      }
    } finally {
      Diagnostic.removeListener(ev3);
    }
    FileHelper.deleteAll(file);
  }
}

