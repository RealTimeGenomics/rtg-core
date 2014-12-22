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

import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.util.zip.CRC32;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.bytecompression.CompressedByteArray;
import com.rtg.util.io.FileUtils;

/**
 * File pair implementation for writing compressed data
 */
@TestClass("com.rtg.reader.SequencesWriterTest")
class CompressedSequenceFilePair extends NormalSequenceFilePair {

  /**
   * Constructor
   * @param dir directory to write files to
   * @param fileNum file number to write
   * @param quality whether to write quality data
   * @param limit largest number of values allowed
   * @param seqValues range of possible values
   * @param checksumSeq checksum tracker for sequence data
   * @param checksumQual checksum tracker for quality data
   * @throws IOException if an IO error occurs
   */
  public CompressedSequenceFilePair(File dir, int fileNum, boolean quality, final long limit, final int seqValues, CRC32 checksumSeq, CRC32 checksumQual) throws IOException {
    super(new FileBitwiseOutputStream(SdfFileUtils.sequenceDataFile(dir, fileNum), CompressedByteArray.minBits(seqValues)),
          quality ? new FileCompressedOutputStream(SdfFileUtils.qualityDataFile(dir, fileNum), SdfWriter.MAX_QUALITY_VALUE) : null,
          new DataOutputStream(FileUtils.createOutputStream(SdfFileUtils.sequencePointerFile(dir, fileNum), false)),
          quality, limit, checksumSeq, checksumQual);

  }
}
