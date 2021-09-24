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
package com.rtg.ngs;


import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.Writer;

import com.rtg.index.hash.ngs.OutputProcessor;
import com.rtg.launcher.HashingRegion;
import com.rtg.reader.NamesInterface;


/**
 * Default implementation for output processing.
 */
public class DefaultOutputProcessor implements OutputProcessor {

  private OutputStreamWriter mOut;
  private final NamesInterface mNames;
  private boolean mHeaderWritten = false;
  private final File mFile;
  private final HashingRegion mRegion;

  /**
   * Construct an output processor
   * @param params parameters
   * @throws IOException if an I/O Error occurs
   */
  public DefaultOutputProcessor(final NgsParams params) throws IOException {
    mOut = new OutputStreamWriter(params.unusedOutStream());
    if (params.searchParams() == null) {
      mNames = null;
    } else {
      mNames = params.searchParams().reader().names();
    }
    mFile = null;
    mRegion = HashingRegion.NONE;
  }
  /**
   * Construct an output processor
   * @param stream output destination
   * @param names name store
   * @param region the region this processor is dealing with
   * @param output the file this processor is outputting to
   */
  public DefaultOutputProcessor(final OutputStream stream, NamesInterface names, HashingRegion region, File output) {
    mOut = new OutputStreamWriter(stream);
    mNames = names;
    mFile = output;
    mRegion = region;
  }

  public File getFile() {
    return mFile;
  }

  public HashingRegion getRegion() {
    return mRegion;
  }

  /**
   * @see OutputProcessor#process(long, java.lang.String, int, int, int, int)
   */
  @Override
  public void process(final long templateId, final String frame, final int readId, final int tStart, final int score, final int scoreIndel) throws IOException {
    final String templateName;
    final boolean writeTemplateName = mNames != null;
    if (!mHeaderWritten) {
      mHeaderWritten = true;
      writeHeader(mOut, writeTemplateName);
    }
    if (writeTemplateName) {
      templateName = mNames.name(templateId);
    } else {
      templateName = String.valueOf(templateId);
    }
    mOut.append("")
      .append(templateName)
      .append("\t")
      .append(frame)
      .append("\t")
      .append(String.valueOf(readId))
      .append("\t")
      .append(String.valueOf(tStart + 1))
      .append("\t")
      .append(String.valueOf(score))
      .append("\t")
      .append(String.valueOf(scoreIndel))
      .append(LS);
  }

  @Override
  public void finish() throws IOException {
    mOut.flush();
  }

  @Override
  @SuppressWarnings("try")
  public void close() throws IOException {
    if (mOut != null) {
      try (final OutputStreamWriter ignored = mOut) {
        mOut = null;
      }
    }
  }


 /**
  * Write the header for result file indicating column names.
  * @param writer destination for header
  * @param templateName Indicates if we are writing Template names or id.
  * @throws IOException if problems when writing output.
  */
  static void writeHeader(final Writer writer, final boolean templateName) throws IOException {
    writer.append("#")
      .append(templateName ? "template-name" : "template-id")
      .append("\tframe\tread-id\ttemplate-start\tscore\tscore-indel").append(LS);
  }

  @Override
  public OutputProcessor threadClone(HashingRegion region) {
    throw new UnsupportedOperationException("DefaultOutputProcessor is not thread-safe");
  }

  @Override
  public void threadFinish() throws IOException {
    try {
      finish();
    } finally {
      close();
    }
  }
}

