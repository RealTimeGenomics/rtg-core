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
public final class SdfWriterWrapper implements WriterWrapper {

  private final boolean mIsPaired;
  private final boolean mHasQuality;
  private final boolean mHasNames;
  private SdfReaderWrapper mReader;
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
    mReader = reader;
    mIsPaired = reader.isPaired();
    mHasQuality = reader.hasQualityData();
    mHasNames = reader.hasNames();

    if (mIsPaired) {
      FileUtils.ensureOutputDirectory(baseDir);
      mSingle = null;
      mLeft = new SdfWriter(new File(baseDir, "left"), Constants.MAX_FILE_SIZE, reader.getPrereadType(), mHasQuality, mHasNames, forceCompression || reader.left().compressed(), reader.type());
      mLeft.setPrereadArm(reader.left().getArm());
      mLeft.setCommandLine(CommandLine.getCommandLine());
      mRight = new SdfWriter(new File(baseDir, "right"), Constants.MAX_FILE_SIZE, reader.getPrereadType(), mHasQuality, mHasNames, forceCompression || reader.right().compressed(), reader.type());
      mRight.setPrereadArm(reader.right().getArm());
      mRight.setCommandLine(CommandLine.getCommandLine());
      mRight.setSdfId(mLeft.getSdfId());
    } else {
      mLeft = null;
      mRight = null;
      mSingle = new SdfWriter(baseDir, Constants.MAX_FILE_SIZE, reader.getPrereadType(), mHasQuality, mHasNames, forceCompression || reader.single().compressed(), reader.type());
      mSingle.setPrereadArm(reader.single().getArm());
      mSingle.setCommandLine(CommandLine.getCommandLine());
    }
    setSdfId(new SdfId());
  }

  /**
   * Set a new reader from which sequences will be selected -- only
   * needed for <code>sdfsplit</code> during merge operation!
   * @param reader the new reader to supply sequence data.
   */
  void setReader(SdfReaderWrapper reader) {
    mReader = reader;
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

  @Override
  public void writeSequence(long seqId, byte[] dataBuffer, byte[] qualityBuffer) throws IllegalStateException, IOException {
    if (mIsPaired) {
      writeSequence(mReader.left(), seqId, mLeft, dataBuffer, qualityBuffer);
      writeSequence(mReader.right(), seqId, mRight, dataBuffer, qualityBuffer);
    } else {
      writeSequence(mReader.single(), seqId, mSingle, dataBuffer, qualityBuffer);
    }
  }

  private void writeSequence(SequencesReader reader, long seqId, SdfWriter writer, byte[] dataBuffer, byte[] qualityBuffer) throws IllegalArgumentException, IllegalStateException, IOException {
    final int length = reader.read(seqId, dataBuffer);
    if (mHasQuality) {
      reader.readQuality(seqId, qualityBuffer);
    }
    writer.startSequence(mHasNames ? reader.fullName(seqId) : null);
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

