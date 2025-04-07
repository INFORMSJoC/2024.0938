#include <sstream>
#include <fstream>
#include <iomanip>
#include "CGStruct.h"

void CGResult::print() const {

    std::cout << "\n" << "Number of points: " << this->n << "\n";
    std::cout << "Number of clusters: " << this->k << "\n";
    std::cout << "Number of components: " << this->n_components << "\n\n";

    std::cout << "Iter: " << this->n_iter << "\n";
    std::cout << "Initial upper bound: " << this->start_ub << "\n";
    std::cout << "Final lower bound: " << this->best_lb << "\n";
    std::cout << "Final upper bound: " << this->best_ub << "\n";
    std::cout << "Initial IP gap: " << (this->start_ub - this->best_lb) / this->start_ub * 100 << " [%]\n";
    std::cout << "Final IP gap: " << (this->best_ub - this->best_lb) / this->best_ub * 100 << " [%]\n\n";

    std::cout << "Initial number of covering constraints: " << this->n_covering_start << "\n";
    std::cout << "Final number of covering constraints: " << this->n_covering_end << "\n";
    std::cout << "Number of times the partition has been updated: " << this->n_update << "\n";
    std::cout << "Count active bounds on lambda: " << this->count_times_active << "\n\n";

    std::cout << "Total time master (sec) = " << this->time_master << "\n";
    std::cout << "Total time pricing (sec) = " << this->time_pricing << "\n";
    std::cout << "Total time (sec) = " <<  this->time << "\n\n";

    std::cout << "Average time master (sec): " << this->avg_time_master << "\n";
    std::cout << "Average time pricing (sec): " << this->avg_time_pricing << "\n";
    std::cout << "Average incompatibilities (best d-incompatible column) = " << this->avg_target_degree << "\n";
    std::cout << "Average incompatibilities (best incompatible column) = " << this->avg_best_degree << "\n";
    std::cout << "Average covering = " << this->avg_covering << "\n";
    std::cout << "Average number of generated columns = " << this->avg_n_incompatible << "\n";
    std::cout << "Average number of compatible columns = " << this->avg_n_compatible << "\n\n";

}

std::string CGResult::to_csv_header() {
    std::stringstream ss;
    ss << "n" << "," << "k" << "," << "n compo" << "," << "UB start" << ",";
    ss << "UB" << "," << "iter" << "," << "LB" << "," << "IP gap [%]" << ",";
    ss << "n start" << "," << "n end" << "," << "avg n" << "," << "update" << "," << "act lambda" << ",";
    ss << "avg inc" << "," << "avg d-inc" << "," << "avg col added" << ",";
    ss << "avg time RMP [s]" << "," << "avg time PP [s]" << ",";
    ss << "time RMP [s]" << "," << "time PP [s]" << "," << "time [s]";
    return ss.str();
}

std::string CGResult::to_csv_line() const {
    std::stringstream ss;
    ss << this->n << "," << this->k << "," << this->n_components << "," << this->start_ub << ",";
    ss << this->best_ub << "," << this->n_iter << "," << this->best_lb << ",";
    ss << std::fixed << std::setprecision(2) << std::abs(this->best_ub - this->best_lb) / this->best_ub * 100 << ",";
    ss << this->n_covering_start << "," << this->n_covering_end << "," << std::fixed << std::setprecision(2) << this->avg_covering << ",";
    ss << this->n_update << "," << this->count_times_active << ",";
    ss << std::fixed << std::setprecision(2) << this->avg_best_degree << ",";
    ss << std::fixed << std::setprecision(2) << this->avg_target_degree << ",";
    ss << std::fixed << std::setprecision(2) << this->avg_n_compatible << ",";
    ss << std::fixed << std::setprecision(2) << this->avg_time_master << ",";
    ss << std::fixed << std::setprecision(2) << this->avg_time_pricing << ",";
    ss << std::fixed << std::setprecision(2) << this->time_master << ",";
    ss << std::fixed << std::setprecision(2) << this->time_pricing << ",";
    ss << std::fixed << std::setprecision(2) << this->time;
    return ss.str();
}

void CGResult::save_assignment(const char *filename) const {
    if (!this->cr_end->assignment.empty()) {
        std::ofstream outFile(filename);
        // the important part
        for (const auto &e : this->cr_end->assignment) outFile << e << "\n";
        outFile.close();
    }
}

CGResult::~CGResult() {

    delete (this->grb_end);
    delete (this->cr_end);
    delete (this->model_end);

}

void CGLog::header() {

    std::cout << "\n" << "    ------------------------------------------------------------------------"
                 "  COLUMN GENERATION WITH DCA FOR THE MSSC PROBLEM  "
                 "--------------------------------------------------------------------\n" <<
              std::setw(8) << "Iter" << " " <<
              std::setw(10) << "m" << " " <<
              std::setw(15) << "Col ARMP" << " " <<
              std::setw(15) << "Opt ARMP" << " " <<
              std::setw(12) << "Time ARMP" << " " <<
              //std::setw(9) << "Inc" << " " <<
              std::setw(8) << "Comp" << " " <<
              std::setw(15) << "Best Inc" << " " <<
              std::setw(15) << "Best Comp" << " " <<
              std::setw(11) << "Ratio" << " " <<
              std::setw(12) << "Time PP" << " " <<
              std::setw(15) << "LB Inc" << " " <<
              //std::setw(11) << "CG Gap [%]" << " " <<
              std::setw(15) << "CG Gap [%]" <<
              std::setw(15) << "IP Gap [%]" <<
              std::setw(15) << "Time [s]" << std::endl;

}

void CGLog::footer() {

    std::cout << "    ----------------------------------------------------------------"
                 "----------------------------------------------------"
                 "---------------------------------------------------------------------------\n\n";

}

void CGLog::print() const {

    if (this->verbose) {
        std::cout <<
                  std::setw(8) << this->iter << " " <<
                  std::setw(10) << this->n_covering << " " <<
                  std::setw(15) << this->col_armp << " " <<
                  std::setw(15) << this->opt_armp << " " <<
                  std::setw(12) << std::fixed << std::setprecision(2) << this->time_armp << " " <<
                  //std::setw(9) << this->n_incompatible << " " <<
                  std::setw(8) << this->n_compatible << " " <<
                  std::setw(15) << this->pricing_incompatible << " " <<
                  std::setw(15) << this->pricing_compatible << " " <<
                  std::setw(11) << std::fixed << std::setprecision(2) << this->col_ratio << " " <<
                  std::setw(12) << std::fixed << std::setprecision(2) << this->time_pricing << " " <<
                  std::setw(15) << this->lb_incomp << " " <<
                  //std::setw(11) << std::fixed << std::setprecision(2) << this->gap * 100 << " " <<
                  std::setw(15) << std::fixed << std::setprecision(3) << this->best_gap * 100 << " " <<
                  std::setw(15) << std::fixed << std::setprecision(3) << this->ip_gap * 100 << " " <<
                  std::setw(15) << std::fixed << std::setprecision(2) << this->time_iter;
        std::cout << "\n";
    }

}
