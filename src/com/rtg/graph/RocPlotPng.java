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
package com.rtg.graph;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import com.reeltwo.plot.Graph2D;
import com.reeltwo.plot.renderer.GraphicsRenderer;
import com.reeltwo.plot.ui.ImageWriter;
import com.rtg.util.io.FileUtils;

/**
 */
public final class RocPlotPng {

  private RocPlotPng() {

  }

  static void rocPngImage(List<File> fileList, List<String> nameList, String title, boolean scores, int lineWidth, File pngFile) throws IOException {
    final Map<String, DataBundle> data = new LinkedHashMap<>();
    for (int i = 0; i < fileList.size(); i++) {
      final DataBundle db = ParseRocFile.loadStream(new NullProgressDelegate(), FileUtils.createInputStream(fileList.get(i), false), nameList.get(i), false);
      data.put(db.getTitle(), db);
    }
    final ImageWriter iw = new ImageWriter(new GraphicsRenderer());

    final ArrayList<String> paths = new ArrayList<>(data.keySet());

    final Graph2D graph = new RocPlot.RocGraph2D(paths, lineWidth, scores, data, title != null ? title : "ROC");
    iw.toPNG(pngFile, graph, 800, 600, null);
  }

  private static class NullProgressDelegate implements ProgressDelegate {
    @Override
    public void setProgress(int progress) {
    }

    @Override
    public void addFile(int numberLines) {
    }

    @Override
    public void done() {
    }
  }
}
