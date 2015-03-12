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
package com.rtg.visualization;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;

import com.rtg.mode.DNA;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.util.diagnostic.NoTalkbackSlimException;

/**
 * Class containing helper utilities for <code>tview</code> and
 * <code>Aview</code>
 *
 */
public final class ReferenceHelper {

  private ReferenceHelper() {
  }

  static final char[] BASES = {Character.toLowerCase(DNA.values()[0].toString().charAt(0)),
    Character.toLowerCase(DNA.values()[1].toString().charAt(0)),
    Character.toLowerCase(DNA.values()[2].toString().charAt(0)),
    Character.toLowerCase(DNA.values()[3].toString().charAt(0)),
    Character.toLowerCase(DNA.values()[4].toString().charAt(0))
  };

  /**
   * Get the reference bases in the given hash map as a string
   *
   * @param template hash map containing template
   * @param ref reference name
   * @param start start inclusive 0 based
   * @param end end of template exclusive
   * @return string containing reference within given range
   */
  public static String getReference(HashMap<String, byte[]> template, String ref, int start, int end) {
    final byte[] b = template.get(ref);
    return getReferenceFromBytes(b, start, end);
  }

  /**
   * Get the reference from bytes
   *
   * @param start start inclusive 0 based
   * @param end end of template exclusive
   * @param b bytes
   * @return string containing reference
   */
  public static String getReferenceFromBytes(byte[] b, int start, int end) {
    final StringBuilder sb = new StringBuilder();
    for (int i = start; i < end && i < b.length; i++) {
      sb.append(ReferenceHelper.BASES[b[i]]);
    }
    return sb.toString();
  }

  static String referenceStringFromBytes(byte[] b, int inclusiveStart, int inclusiveEnd) {
    final StringBuilder sb = new StringBuilder();
    for (int i = inclusiveStart; i <= inclusiveEnd && i < b.length; i++) {
      sb.append(ReferenceHelper.BASES[b[i]]);
    }
    return sb.toString();
  }

  /**
   * Load a template in given hash map
   *
   * @param templateFile SDF containing template
   * @param template hash map
   * @param err error to print statistics on, can be <code>null</code>
   * @throws IOException if error
   */
  public static void loadTemplate(final File templateFile, final HashMap<String, byte[]> template, final PrintStream err) throws IOException {
    try (SequencesReader reader = SequencesReaderFactory.createDefaultSequencesReaderCheckEmpty(templateFile)) {
      for (long seq = 0; seq < reader.numberSequences(); seq++) {
        template.put(reader.name(seq), reader.read(seq));
      }
      if (err != null) {
        err.println("Template loaded : " + reader.numberSequences());
      }
    }
  }

  /**
   * Loads a template in hash map
   *
   * @param dir SDF containing template
   * @return hash map containing template
   * @throws IOException if error
   */
  public static HashMap<String, byte[]> loadAllTemplate(File dir) throws IOException {
    final HashMap<String, byte[]> template = new HashMap<>();
    loadTemplate(dir, template, null);
    return template;
  }

  static byte[] loadSingleTemplate(File templateFile, String sequenceName) throws IOException {
    try (SequencesReader reader = SequencesReaderFactory.createDefaultSequencesReaderCheckEmpty(templateFile)) {
      for (long seq = 0; seq < reader.numberSequences(); seq++) {
        if (reader.name(seq).equals(sequenceName)) {
          return reader.read(seq);
        }
      }
    }
    throw new NoTalkbackSlimException("Given sequence name not found : " + sequenceName);
  }

  static int templateLength(File templateFile, String sequenceName) throws IOException {
    try (SequencesReader reader = SequencesReaderFactory.createDefaultSequencesReaderCheckEmpty(templateFile)) {
      for (long seq = 0; seq < reader.numberSequences(); seq++) {
        if (reader.name(seq).equals(sequenceName)) {
          return reader.length(seq);
        }
      }
      throw new NoTalkbackSlimException("Given sequence name not found : " + sequenceName);
    }
  }

  static byte[] loadReference(File templateFile, String sequenceName, int start, int len) throws IOException {
    try (SequencesReader reader = SequencesReaderFactory.createDefaultSequencesReaderCheckEmpty(templateFile)) {
      for (long seq = 0; seq < reader.numberSequences(); seq++) {
        if (reader.name(seq).equals(sequenceName)) {
          final byte[] b = new byte[len]; // XXX What happens if start + len > length(seq) ?
          reader.read(seq, b, start, len);
          return b;
        }
      }
      throw new NoTalkbackSlimException("Given sequence name not found : " + sequenceName);
    }
  }
}
