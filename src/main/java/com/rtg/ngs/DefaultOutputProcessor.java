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

