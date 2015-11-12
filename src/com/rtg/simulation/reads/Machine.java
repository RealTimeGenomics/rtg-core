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
package com.rtg.simulation.reads;

import java.io.IOException;

import com.rtg.reader.PrereadType;
import com.rtg.reader.SdfId;
import com.rtg.util.machine.MachineType;

/**
 * Interface to represent different machines for read simulation.
 */
public interface Machine {

  /**
   * Sets the minimum and maximum quality values that will be output for bases.
   *
   * @param minq the minimum quality value permitted.
   * @param maxq the maximum quality value permitted.
   */
  void setQualRange(byte minq, byte maxq);

  /**
   * Sets the destination <code>ReadWriter</code> to which the simulated reads will be sent.
   *
   * @param rw a <code>ReadWriter</code> value
   */
  void setReadWriter(ReadWriter rw);

  /**
   * Specifies the list of template sets used during generation.
   * @param templateIds an array containing an ID for each template set
   */
  void identifyTemplateSet(SdfId... templateIds);

  /**
   * Specifies the original reference template used during mutated genome generation.
   * @param referenceId the ID of the original reference template.
   */
  void identifyOriginalReference(SdfId referenceId);

  /**
   * Take a fragment and generate a read.
   * @param id fragment name
   * @param fragmentStart 0 based start position of fragment within sequence. Negative values are valid, but shouldn't be used for any calculations within this method..
   * @param data residues
   * @param length amount of data that is valid
   * @throws IOException guess
   */
  void processFragment(String id, int fragmentStart, byte[] data, int length) throws IOException;

  /**
   * Total residues emitted in generated reads
   * @return the number
   */
  long residues();


  /**
   * @return true if this machine produces paired end data.
   */
  boolean isPaired();

  /**
   * @return type to use when writing generated reads to SDF
   */
  PrereadType prereadType();

  /**
   * @return which type of machine is being simulated.
   */
  MachineType machineType();

  /**
   * @return a textual representation summary of the actions histogram
   */
  String formatActionsHistogram();
}
