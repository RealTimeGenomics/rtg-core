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
import com.rtg.mode.SequenceType;
import com.rtg.util.io.FileUtils;

/**
 * Base class to make implementation of readers easier
 */
@TestClass(value = {"com.rtg.reader.DefaultSequencesReaderTest", "com.rtg.reader.CompressedMemorySequencesReaderTest"})
public abstract class AbstractSequencesReader implements AnnotatedSequencesReader {

  private static final String README_FILENAME = "readme.txt";

  protected abstract IndexFile index();

  @Override
  public String currentFullName() throws IllegalStateException, IOException {
    return currentName() + currentNameSuffix();
  }

  @Override
  public String fullName(long sequenceIndex) throws IOException {
    return name(sequenceIndex) + nameSuffix(sequenceIndex);
  }

  @Override
  public long[] residueCounts() {
    return index().getResidueCounts();
  }

  @Override
  public long dataChecksum() {
    return index().getDataChecksum();
  }
  @Override
  public long qualityChecksum() {
    return index().getQualityChecksum();
  }
  @Override
  public long nameChecksum() {
    return index().getNameChecksum();
  }


  @Override
  public long[] histogram() {
    return index().getNHistogram();
  }

  @Override
  public long[] posHistogram() {
    return index().getPosHistogram();
  }

  @Override
  public double globalQualityAverage() {
    return index().getQSAverage();
  }

  @Override
  public boolean hasHistogram() {
    return index().hasNStats();
  }

  @Override
  public long longestNBlock() {
    return index().getLongestNBlock();
  }

  @Override
  public long nBlockCount() {
    return index().getNBlockCount();
  }

  @Override
  public PrereadArm getArm() {
    return index().getPrereadArm();
  }

  @Override
  public PrereadType getPrereadType() {
    return index().getPrereadType();
  }

  @Override
  public SdfId getSdfId() {
    return index().getSdfId();
  }

  @Override
  public double[] positionQualityAverage() {
    return index().getQSPositionAverageHistogram();
  }

  @Override
  public long sdfVersion() {
    return index().getVersion();
  }

  @Override
  public String comment() {
    return index().getComment();
  }

  @Override
  public String commandLine() {
    return index().getCommandLine();
  }

  @Override
  public String samReadGroup() {
    return index().getSamReadGroup();
  }

  @Override
  public boolean compressed() {
    return index().getSequenceEncoding() == IndexFile.SEQUENCE_ENCODING_COMPRESSED;
  }

  @Override
  public long suffixChecksum() {
    return index().getNameSuffixChecksum();
  }

  @Override
  public long totalLength() {
    return index().getTotalLength();
  }

  @Override
  public long maxLength() {
    return index().getMaxLength();
  }

  @Override
  public long minLength() {
    return index().getMinLength();
  }

  @Override
  public boolean hasNames() {
    return index().hasNames();
  }

  @Override
  public boolean hasQualityData() {
    return index().hasQuality();
  }

  @Override
  public SequenceType type() {
    return SequenceType.values()[index().getSequenceType()];
  }

  @Override
  public String getReadMe() throws IOException {
    if (path() == null) {
      return null;
    }
    final File readMe = new File(path(), README_FILENAME);
    if (!readMe.isFile()) {
      return null;
    }
    return FileUtils.fileToString(readMe);
  }
}
