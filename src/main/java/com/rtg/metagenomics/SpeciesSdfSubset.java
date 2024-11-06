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
package com.rtg.metagenomics;

import static com.rtg.launcher.CommonFlags.INPUT_FLAG;
import static com.rtg.launcher.CommonFlags.OUTPUT_FLAG;
import static com.rtg.util.cli.CommonFlagCategories.FILTERING;
import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Map;

import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.LoggedCli;
import com.rtg.reader.AbstractSdfWriter.SequenceNameHandler;
import com.rtg.reader.Label;
import com.rtg.reader.ReaderUtils;
import com.rtg.reader.SdfReaderWrapper;
import com.rtg.reader.SdfWriterWrapper;
import com.rtg.taxonomy.Taxonomy;
import com.rtg.taxonomy.TaxonomyUtils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.diagnostic.Warnings;
import com.rtg.util.io.LogStream;

/**
 * Pulls out a subset of sequences from one SDF into another SDF.
 *
 */
public final class SpeciesSdfSubset extends LoggedCli {

  private static final String TAXONOMY_FILE_FLAG = "taxonomy-file";

  private static final Validator VALIDATOR = new Validator() {
    @Override
    public boolean isValid(CFlags flags) {
      return CommonFlags.validateOutputDirectory(flags);
    }
  };

  private byte[] mData = null;
  private byte[] mQualities = null;
  private final Warnings mWarnings = new Warnings();
  private final Warnings.Warning mInvalidId = mWarnings.create(5, "Invalid sequence name", false);
  private long mWritten = 0;
  private SdfWriterWrapper mWriter = null;
  private SdfReaderWrapper mReader = null;
  private Map<String, Long> mNames = null;

  @Override
  public String moduleName() {
    return "speciessdfsubset";
  }

  @Override
  public String description() {
    return null;
  }

  @Override
  protected void initFlags() {
    CommonFlagCategories.setCategories(mFlags);
    mFlags.setDescription("Extracts a subset of sequences from one SDF and outputs them to another SDF.");
    mFlags.registerRequired('i', INPUT_FLAG, File.class, CommonFlags.SDF, "input SDF").setCategory(INPUT_OUTPUT);
    mFlags.registerRequired('o', OUTPUT_FLAG, File.class, CommonFlags.SDF, "output SDF").setCategory(INPUT_OUTPUT);
    mFlags.registerRequired('t', TAXONOMY_FILE_FLAG, File.class, CommonFlags.FILE, "file containing taxonomy").setCategory(FILTERING);
    mFlags.setValidator(VALIDATOR);
  }

  @Override
  protected File outputDirectory() {
    return (File) mFlags.getValue(OUTPUT_FLAG);
  }

  private void warnInvalidId(String seqid) {
    mInvalidId.warn(seqid);
  }

  private void getSequence(long seqid) throws IOException {
    if ((seqid < 0) || (seqid >= mReader.numberSequences())) {
      warnInvalidId(String.valueOf(seqid));
      return;
    }
    final int length = mReader.maxLength();
    if (mData == null || mData.length < length) {
      mData = new byte[length];
      if (mReader.hasQualityData()) {
        mQualities = new byte[length];
      }
    }
    mWriter.writeSequence(seqid, mData, mQualities);

    if (++mWritten % 1000 == 0) {
      Diagnostic.progress("Extracted " + mWritten + " sequences");
    }
  }

  private void getSequence(String seqName) throws IOException {
    final Long id = mNames.get(seqName);
    if (id != null) {
      getSequence(id);
    } else {
      warnInvalidId(seqName);
    }
  }

  @Override
  protected int mainExec(OutputStream out, LogStream log) throws IOException {
    final File inputSdf = (File) mFlags.getValue(INPUT_FLAG);
    mReader = new SdfReaderWrapper(inputSdf, false, false);
    if (!mReader.hasNames()) {
      throw new NoTalkbackSlimException("Names are required in SDF.");
    }

    mNames = ReaderUtils.getSequenceNameMap(mReader.isPaired() ? mReader.left() : mReader.single());

    final File output = (File) mFlags.getValue(OUTPUT_FLAG);
    mWriter = new SdfWriterWrapper(output, mReader, true);

    final Taxonomy taxonomy = new Taxonomy();
    final File taxFile = (File) mFlags.getValue(TAXONOMY_FILE_FLAG);
    try (FileInputStream fis = new FileInputStream(taxFile)) {
      taxonomy.read(fis);
    }

    final SequenceNameHandler snh = new SequenceNameHandler();

    final File taxonomyToNameFile = new File(inputSdf, TaxonomyUtils.TAXONOMY_TO_SEQUENCE_FILE);
    final ArrayList<String> taxIdsToSeqNames = new ArrayList<>();
    try (BufferedReader br = new BufferedReader(new FileReader(taxonomyToNameFile))) {
      String line;
      while ((line = br.readLine()) != null) {
        line = line.trim();
        if (line.startsWith("#") || line.length() == 0) {
          taxIdsToSeqNames.add(line);
          continue;
        }
        if (line.length() > 0) {
          final String[] parts = line.split("\t");
          if (parts.length != 2) {
            throw new IOException("File:" + taxonomyToNameFile.getCanonicalPath() + " should be two fields in line: " + line);
          }
          final int taxId = Integer.parseInt(parts[0]);
          if (taxonomy.contains(taxId)) {
            taxIdsToSeqNames.add(line);
            final Label label = snh.handleSequenceName(parts[1]);
            getSequence(label.label());
          }
        }
      }
    }
    mWriter.close();
    mReader.close();
    Diagnostic.progress("Extracted " + mWritten + " sequences");

    final File lookupFile = new File(output, TaxonomyUtils.TAXONOMY_TO_SEQUENCE_FILE);
    Diagnostic.progress("Started writing " + lookupFile.getPath());
    try (BufferedWriter bw = new BufferedWriter(new FileWriter(lookupFile))) {
      for (final String line : taxIdsToSeqNames) {
        bw.write(line);
        bw.newLine();
      }
    }
    Diagnostic.progress("Finished writing " + lookupFile.getPath());

    final File taxonomyFile = new File(output, TaxonomyUtils.TAXONOMY_FILE);
    Diagnostic.progress("Started writing " + taxonomyFile.getPath());
    try (BufferedWriter bw = new BufferedWriter(new FileWriter(taxonomyFile))) {
      taxonomy.write(bw);
    }
    Diagnostic.progress("Finished writing " + taxonomyFile.getPath());
    mWarnings.report();

    return 0;
  }

  /**
   * Main method for pulling out reads into another SDF.
   * @param args arguments for <code>SdfSubset</code>
   */
  public static void main(String[] args) {
    new SpeciesSdfSubset().mainExit(args);
  }
}
