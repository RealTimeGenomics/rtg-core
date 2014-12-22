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
package com.rtg.tabix;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PipedInputStream;
import java.io.PipedOutputStream;

import com.rtg.sam.BamIndexer;
import com.rtg.util.IORunnable;
import com.rtg.util.IORunnableProxy;
import com.rtg.util.io.FileUtils;

/**
 * Gets an output stream that has a parallel thread performing tabix/BAM indexing
 */
public final class IndexingStreamCreator implements Closeable {

  private final File mOutputFile;
  private final boolean mCompress;
  private final boolean mCreateIndexIfPossible;
  private final TabixIndexer.IndexerFactory mIndexerFactory;

  private IORunnableProxy mProxy;
  private Thread mIndexThread;
  private OutputStream mOutputStream;

  /**
   * Construct the stream creator
   * @param outputFile the output file to be written to
   * @param stdout the output stream of stdout (will only be used if output file is null)
   * @param compress true if the output should be gzip compressed (also true if bam is true)
   * @param indexerFactory the indexer factory corresponding to the format being written. Set this to null to indicate BAM indexing.
   * @param createIndexIfPossible true if an index should be created
   */
  public IndexingStreamCreator(File outputFile, OutputStream stdout, boolean compress, TabixIndexer.IndexerFactory indexerFactory, boolean createIndexIfPossible) {
    if (indexerFactory == null) {
      if (!compress) {
        throw new IllegalArgumentException("Uncompressed BAM not supported");
      }
      if (outputFile == null) {
        throw new IllegalArgumentException("BAM output to stdout not supported");
      }
    }
    mOutputFile = outputFile;
    mOutputStream = mOutputFile != null ? null : stdout;
    mCompress = compress;
    mIndexerFactory = indexerFactory;
    mCreateIndexIfPossible = createIndexIfPossible;
  }

  /**
   * @return true if output will go to a file, false if output will go to stdout.
   */
  public boolean isFileDestination() {
    return mOutputFile != null;
  }

  @Override
  public void close() throws IOException {
    try {
      if (mOutputStream != null) {
        mOutputStream.close();
      }
    } finally {
      if (mIndexThread != null) {
        try {
          mIndexThread.join();
        } catch (InterruptedException e) {
          throw new IOException("Execution was interrupted", e);
        } finally {
          mProxy.checkError();
        }
      }
    }
  }

  /**
   * Create the output stream.
   * @return the output stream
   * @throws IOException if there is a problem
   */
  public OutputStream createStreamsAndStartThreads() throws IOException {
    if (mIndexerFactory == null) {
      throw new UnsupportedOperationException("This method is not for use with BAM writing");
    }
    return createStreamsAndStartThreads(-1, true, false);
  }

  /**
   * Create the output stream
   * @param numberReferences the number of reference sequences (needed for BAM)
   * @param expectHeader true if BAM indexing should expect to see a header at the start of the file (may be false if doing indexing of chunks)
   * @param terminateBlockGzip true if the output stream should contain a termination block (may be false if doing indexing of chunks)
   * @return the output stream
   * @throws IOException if there is a problem
   */
  public OutputStream createStreamsAndStartThreads(int numberReferences, boolean expectHeader, boolean terminateBlockGzip) throws IOException {
    if (mOutputStream == null) {
      final boolean bam = mIndexerFactory == null;
      final boolean gzonly = !bam && mCompress;
      if (mCreateIndexIfPossible && mCompress) {
        final PipedInputStream pipeToIndexIn = new PipedInputStream(); //closed by IndexRunner
        final PipedOutputStream pipeToIndexOut = new PipedOutputStream(pipeToIndexIn);
        mOutputStream = FileUtils.createTeedOutputStream(mOutputFile, pipeToIndexOut, gzonly, false, terminateBlockGzip);
        final OutputStream indexOutputStream;
        if (bam) {
          indexOutputStream = FileUtils.createOutputStream(BamIndexer.indexFileName(mOutputFile));
        } else {
          assert mCompress;
          indexOutputStream = FileUtils.createOutputStream(TabixIndexer.indexFileName(mOutputFile));
        }
        final IndexRunner indexRunner = new IndexRunner(pipeToIndexIn, indexOutputStream, mIndexerFactory, expectHeader, numberReferences, mOutputFile.toString());
        mProxy = new IORunnableProxy(indexRunner);
        mIndexThread = new Thread(mProxy);
        mIndexThread.start();
      } else {
        mOutputStream = FileUtils.createOutputStream(mOutputFile, gzonly, false, true);
      }
    }
    return mOutputStream;
  }

  /**
   * Runner for producing the tabix index or BAM file index from an input stream.
   */
  public static class IndexRunner implements IORunnable {

    private final String mFilename;
    private final InputStream mStreamToIndex;
    private final OutputStream mIndexOut;
    private final boolean mBam;
    private final boolean mExpectHeader;
    private final int mNumReferences;
    private final TabixIndexer.IndexerFactory mIndexerFactory;

    /**
     * Constructor for index runner.
     * @param streamToIndex the input stream of the file to be indexed.
     * @param indexOut the output stream for the index file.
     * @param indexerFactory the indexer factory corresponding to the output format. Set this to null to indicate BAM indexing.
     * @param expectHeader true if BAM indexing is expecting the input file to have a header, false otherwise
     * @param numReferences the number of reference sequences, used with BAM input with no header
     * @param filename the name of the file being indexed
     */
    public IndexRunner(InputStream streamToIndex, OutputStream indexOut, TabixIndexer.IndexerFactory indexerFactory, boolean expectHeader, int numReferences, String filename) {
      mStreamToIndex = streamToIndex;
      mIndexOut = indexOut;
      mIndexerFactory = indexerFactory;
      mBam = indexerFactory == null;
      mExpectHeader = expectHeader;
      mNumReferences = numReferences;
      mFilename = filename;
    }

    @Override
    public void run() throws IOException {
      try {
        if (mBam) {
          BamIndexer.saveBamIndexNoHeader(mStreamToIndex, mIndexOut, mExpectHeader, mNumReferences);
        } else {
          new TabixIndexer(mStreamToIndex, mIndexOut).saveIndex(mIndexerFactory);
        }
      } catch (final UnindexableDataException e) {
        throw new IOException("Cannot produce index for: " + mFilename + " (try disabling indexing)", e);
      } finally {
        mStreamToIndex.close();
        mIndexOut.close();
      }
    }

  }
}
