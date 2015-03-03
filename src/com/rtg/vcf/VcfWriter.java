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

package com.rtg.vcf;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

import com.rtg.tabix.IndexingStreamCreator;
import com.rtg.tabix.TabixIndexer;
import com.rtg.util.ByteUtils;
import com.rtg.vcf.header.VcfHeader;


/**
 * Writer to write VCF records out into a <code>.vcf</code> output stream.
 *
 */
public class VcfWriter implements Closeable {

  private final IndexingStreamCreator mIndexer;
  private final OutputStream mOut;
  private final VcfHeader mHeader;
  private boolean mHeaderWritten = false;

  /**
   * Creates a new VCF writer, using on-the-fly indexing.
   * @param header header for the file
   * @param outputFile the output file to be written to
   * @param stdout the output stream of stdout (will only be used if output file is null)
   * @param compress true if the output should be gzip compressed (also true if bam is true)
   * @param createIndexIfPossible true if an index should be created
   * @throws java.io.IOException if there is a problem during writing.
   */
  public VcfWriter(VcfHeader header, File outputFile, OutputStream stdout, boolean compress, boolean createIndexIfPossible) throws IOException {
    if (header == null) {
      throw new NullPointerException("header cannot be null");
    }
    mIndexer = new IndexingStreamCreator(outputFile, stdout, compress, new TabixIndexer.VcfIndexerFactory(), createIndexIfPossible);
    mOut = mIndexer.createStreamsAndStartThreads();
    mHeader = header;
  }

  /**
   * create a new VCF writer
   *
   * @param header header for the file
   * @param out stream to write to
   */
  public VcfWriter(VcfHeader header, OutputStream out) {
    if (out == null) {
      throw new NullPointerException("output stream cannot be null");
    }
    if (header == null) {
      throw new NullPointerException("header cannot be null");
    }
    mIndexer = null;
    mOut = out;
    mHeader = header;
  }

  /**
   * @return current header
   */
  public VcfHeader getHeader() {
    return mHeader;
  }

  /**
   * write current header to output stream
   * @throws java.io.IOException if there is an I/O problem
   */
  protected void writeHeader() throws IOException {
    mOut.write(mHeader.toString().getBytes());
  }

  /**
   * Write record
   *
   * @param record record to write
   * @throws IOException if error
   */
  public void write(VcfRecord record) throws IOException {
    if (!mHeaderWritten) {
      mHeaderWritten = true;
      writeHeader();
    }
    writeToStream(record);
  }

  private void writeToStream(VcfRecord record) throws IOException {
    //System.err.println(record.toString());
    mOut.write(record.toString().getBytes());
    ByteUtils.writeNewline(mOut);
  }

  @Override
  public void close() throws IOException {
    if (!mHeaderWritten) {
      mHeaderWritten = true;
      writeHeader();
    }
    mOut.close();
    if (mIndexer != null) {
      mIndexer.close();
    }
  }

}
