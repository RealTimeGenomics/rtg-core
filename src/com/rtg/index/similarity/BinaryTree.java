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
package com.rtg.index.similarity;

import java.io.IOException;
import java.io.PrintStream;
import java.text.DecimalFormat;

import com.rtg.util.StringUtils;

/**
 * A binary tree holding the results of joining.
 *
 */
public class BinaryTree {

  static final DecimalFormat NF = new DecimalFormat("0.0000");

  final BinaryTree mLeft;
  final BinaryTree mRight;
  final double mLeftDistance;
  final double mRightDistance;
  final String mLabel;
  int mNodeNumber;

  BinaryTree(final BinaryTree left, final BinaryTree right, final double leftDistance, final double rightDistance, final String label) {
    mLeft = left;
    mRight = right;
    mLeftDistance = leftDistance;
    mRightDistance = rightDistance;
    mLabel = label;
  }

  double dump(final String indent, final PrintStream ps) {
    double totalDistance = 0;
    if (mLeft != null) {
      totalDistance += mLeft.dump(indent + " ", ps);
    }
    ps.println(indent + mLabel);
    totalDistance += mLeftDistance;
    totalDistance += mRightDistance;
    if (mRight != null) {
      totalDistance += mRight.dump(indent + " ", ps);
    }
    return totalDistance;
  }

  private int nodeNames(final StringBuilder names, final int nodeNumber) {
    int nodeNmb = nodeNumber;
    // In order traversal
    if (mLeft != null) {
      nodeNmb = mLeft.nodeNames(names, nodeNmb);
    }
    mNodeNumber = nodeNmb++;
    names.append("node ")
      .append(mNodeNumber)
      .append(" \"")
      .append(mLabel).append("\"").append(StringUtils.LS);
    if (mRight != null) {
      nodeNmb = mRight.nodeNames(names, nodeNmb);
    }
    return nodeNmb;
  }

  private void edges(final StringBuilder s) {
    // In order traversal
    if (mLeft != null) {
      s.append("edge ").append(mNodeNumber).append("-").append(mLeft.mNodeNumber).append(StringUtils.LS);
      mLeft.edges(s);
    }
    if (mRight != null) {
      s.append("edge ").append(mNodeNumber).append("-").append(mRight.mNodeNumber).append(StringUtils.LS);
      mRight.edges(s);
    }
  }

  private synchronized void toString(final StringBuilder sb) {
    nodeNames(sb, 0);
    edges(sb);
  }

  @Override
  public String toString() {
    final StringBuilder sb = new StringBuilder();
    toString(sb);
    return sb.toString();
  }

  String getLabel() {
    return mLabel;
  }

  private void newick(final int indent, final Appendable out) throws IOException {
    for (int k = 0; k < indent; k++) {
      out.append(' ');
    }
    if (mLeft == null) {
      assert mRight == null;
      out.append(mLabel);
    } else {
      out.append("(").append(StringUtils.LS);
      mLeft.newick(indent + 1, out);
      out.append(':');
      out.append(NF.format(mLeftDistance));
      out.append(",").append(StringUtils.LS);
      mRight.newick(indent + 1, out);
      out.append(':');
      out.append(NF.format(mRightDistance));
      out.append(StringUtils.LS);
      for (int k = 0; k < indent; k++) {
        out.append(' ');
      }
      out.append(")");
    }
  }

  /**
   * Return the tree in New Hampshire format.  Note internal nodes
   * are unlabelled in this format.
   * @param out where to write the tree.
   * @throws IOException if error while writing to out.
   */
  public void newick(final Appendable out) throws IOException {
    newick(0, out);
    out.append(StringUtils.LS);
  }

  private String spaces(final int indent) {
    final StringBuilder sb = new StringBuilder();
    for (int k = 0; k < indent; k++) {
      sb.append(' ');
    }
    return sb.toString();
  }

  private void phyloXml(final int indent, final Appendable out, final double len) throws IOException {
    final String sp = spaces(indent);
    out.append(sp);
    if (Double.isNaN(len)) {
      out.append("<clade>");
    } else {
      out.append("<clade branch_length=\"");
      out.append(NF.format(len));
      out.append("\">");
    }
    out.append(StringUtils.LS);
    if (mLeft == null) {
      assert mRight == null;
      out.append(sp);
      out.append("  <name>");
      out.append(StringUtils.xmlProtect(mLabel));
      out.append("</name>");
      out.append(StringUtils.LS);
    } else {
      mLeft.phyloXml(indent + 2, out, mLeftDistance);
      mRight.phyloXml(indent + 2, out, mRightDistance);
    }
    out.append(sp);
    out.append("</clade>");
    out.append(StringUtils.LS);
  }

  /**
   * Return the tree in XML format.
   * @param out where to write the tree.
   * @throws IOException if error while writing to out.
   */
  public void phyloXml(final Appendable out) throws IOException {
    out.append("<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
    out.append(StringUtils.LS);
    out.append("<phyloxml xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://www.phyloxml.org http://www.phyloxml.org/1.10/phyloxml.xsd\" xmlns=\"http://www.phyloxml.org\">");
    out.append(StringUtils.LS);
    out.append("  <phylogeny rooted=\"false\">");
    out.append(StringUtils.LS);
    out.append("    <description>Tree produced by the RTG phylogeny module</description>");
    out.append(StringUtils.LS);
    phyloXml(4, out, Double.NaN);
    out.append("  </phylogeny>");
    out.append(StringUtils.LS);
    out.append("</phyloxml>");
    out.append(StringUtils.LS);
  }
}
