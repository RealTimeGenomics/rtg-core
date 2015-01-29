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
package com.rtg.reader;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.Constants;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.io.FileUtils;

/**
 * Wrapper for SDF writers that can handle single end or paired end reads
 *
 */
@TestClass(value = {"com.rtg.reader.SdfSplitterTest", "com.rtg.ngs.DummySdfOutputProcessorTest"})
public final class SdfWriterWrapper implements Closeable {

  private final boolean mIsPaired;
  private final boolean mHasQuality;
  final SdfWriter mLeft;
  final SdfWriter mRight;
  final SdfWriter mSingle;

  /**
   * Convenience wrapper for writing.
   * @param baseDir base directory.
   * @param reader the reader that this writer is writing from.
   * @param forceCompression Force SDF to be compressed
   */
  public SdfWriterWrapper(File baseDir, SdfReaderWrapper reader, boolean forceCompression) {
    assert reader != null;
    mIsPaired = reader.isPaired();
    mHasQuality = reader.hasQualityData();
    final boolean hasNames = reader.hasNames();

    if (mIsPaired) {
      FileUtils.ensureOutputDirectory(baseDir);
      mSingle = null;
      mLeft = new SdfWriter(new File(baseDir, "left"), Constants.MAX_FILE_SIZE, reader.getPrereadType(), mHasQuality, hasNames, forceCompression || reader.left().compressed(), reader.type());
      mLeft.setPrereadArm(reader.left().getArm());
      mLeft.setCommandLine(CommandLine.getCommandLine());
      mRight = new SdfWriter(new File(baseDir, "right"), Constants.MAX_FILE_SIZE, reader.getPrereadType(), mHasQuality, hasNames, forceCompression || reader.right().compressed(), reader.type());
      mRight.setPrereadArm(reader.right().getArm());
      mRight.setCommandLine(CommandLine.getCommandLine());
      mRight.setSdfId(mLeft.getSdfId());
    } else {
      mLeft = null;
      mRight = null;
      mSingle = new SdfWriter(baseDir, Constants.MAX_FILE_SIZE, reader.getPrereadType(), mHasQuality, hasNames, forceCompression || reader.single().compressed(), reader.type());
      mSingle.setPrereadArm(reader.single().getArm());
      mSingle.setCommandLine(CommandLine.getCommandLine());
    }
  }

  /**
   * Convenience method.
   * @param id the identifier.
   */
  public void setSdfId(SdfId id) {
    if (mIsPaired) {
      mLeft.setSdfId(id);
      mRight.setSdfId(id);
    } else {
      mSingle.setSdfId(id);
    }
  }

  /**
   * Method to write the current reader sequence to file.
   * @param reader the reader being read from.
   * @param dataBuffer an appropriately sized buffer for data.
   * @param qualityBuffer an appropriately sized buffer for quality.
   * @throws IllegalStateException whenever
   * @throws IOException whenever
   */
  public void writeCurrentSequence(SdfReaderWrapper reader, byte[] dataBuffer, byte[] qualityBuffer) throws IllegalStateException, IOException {
    if (mIsPaired) {
      writeSequence(reader.left(), mLeft, reader.currentFullName(), dataBuffer, qualityBuffer);
      writeSequence(reader.right(), mRight, reader.currentRightFullName(), dataBuffer, qualityBuffer);
    } else {
      writeSequence(reader.single(), mSingle, reader.currentFullName(), dataBuffer, qualityBuffer);
    }
  }

  private void writeSequence(SequencesReader reader, SdfWriter writer, String name, byte[] dataBuffer, byte[] qualityBuffer) throws IllegalArgumentException, IllegalStateException, IOException {
    final int length = reader.readCurrent(dataBuffer);
    if (mHasQuality) {
      reader.readCurrentQuality(qualityBuffer);
    }
    writer.startSequence(name);
    writer.write(dataBuffer, mHasQuality ? qualityBuffer : null, length);
    writer.endSequence();
  }

  /**
   * Function to copy the source templates file into the resulting SDF files
   * if it exists in the original SDF file.
   * @param reader the reader to get the original list from.
   * @throws IOException if an I/O error occurs during reading or writing.
   */
  public void copySourceTemplatesFile(SdfReaderWrapper reader) throws IOException {
    if (mIsPaired) {
      SourceTemplateReadWriter.writeTemplateMappingFile(mLeft.directory(), SourceTemplateReadWriter.readTemplateMap(reader.left().path()));
      SourceTemplateReadWriter.writeTemplateMappingFile(mRight.directory(), SourceTemplateReadWriter.readTemplateMap(reader.right().path()));
      SourceTemplateReadWriter.copyMutationMappingFile(reader.left().path(), mLeft.directory());
      SourceTemplateReadWriter.copyMutationMappingFile(reader.right().path(), mRight.directory());
    } else {
      SourceTemplateReadWriter.writeTemplateMappingFile(mSingle.directory(), SourceTemplateReadWriter.readTemplateMap(reader.single().path()));
      SourceTemplateReadWriter.copyMutationMappingFile(reader.single().path(), mSingle.directory());
    }
  }

  /**
   * Close method for the writers.
   * @throws IOException whenever
   */
  @Override
  public void close() throws IOException {
    if (mLeft != null) {
      mLeft.close();
    }
    if (mRight != null) {
      mRight.close();
    }
    if (mSingle != null) {
      mSingle.close();
    }
  }
}

