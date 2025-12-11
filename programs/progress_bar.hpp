#ifndef PROGRESS_BAR_H_
#define PROGRESS_BAR_H_

void displayProgressBar(int current, int total) {
  int barWidth = 50; // Width of the progress bar
  float progress = (float)current / total;
  int position = barWidth * progress;

  std::cout << "\r[";
  for (int i = 0; i < barWidth; i++) {
    if (i < position) {
      std::cout << "=";
    } else if (i == position) {
      std::cout << ">";
    } else {
      std::cout << " ";
    }
  }
  std::cout << "] " << int(progress * 100.0) << "%";
  std::cout.flush(); // Ensure the output is displayed immediately
}

#endif
