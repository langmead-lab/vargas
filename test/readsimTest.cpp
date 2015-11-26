//
// Created by gaddra on 11/24/15.
//

/**
TEST(readsimTest, test1) {
  vargas::Graph g("data/r5", "data/v5", "out");
  g.exportDOT("out.dot");
  vargas::ReadSim sim(g);
  sim.setReadLen(3);
  sim.addProfile(".*9.*", "out.reads");
  sim.addProfile(".*10.*", "out2.reads");
  sim.setNumReads(10);
  while (sim.updateRead());

  vargas::ReadFile rf("out.reads");
  while (rf.updateRead()) {
//    std::cout << rf.get() << std::endl;
  }
}

 **/