#include "include/utils_export.cuh"

std::string to_string(Experiment::Operation operation) {
    switch (operation) {
        case Experiment::Operation::ADD_SCALAR:
            return "ADD_SCALAR";
        case Experiment::Operation::ADD:
            return "ADD";
        case Experiment::Operation::MUL:
            return "MUL";
        case Experiment::Operation::MUL_SCALAR:
            return "MUL_SCALAR";
        case Experiment::Operation::ROT:
            return "ROT";
        case Experiment::Operation::FAST_ROT:
            return "FAST_ROT";
        case Experiment::Operation::ENCODE:
            return "ENCODE";
        case Experiment::Operation::ENCRYPT:
            return "ENCRYPT";
        case Experiment::Operation::DECRYPT:
            return "DECRYPT";
        case Experiment::Operation::GEN:
            return "GEN";
        case Experiment::Operation::INIT_DC:
            return "INIT_DC";
        case Experiment::Operation::SYNC:
            return "SYNC";
        case Experiment::Operation::PREPROC_DC:
            return "PREPROC_DC";
        case Experiment::Operation::EVAL_DC:
            return "EVAL_DC";
        case Experiment::Operation::EVAL_OP:
            return "EVAL_OP";
        case Experiment::Operation::PRUNING:
            return "PRUNING";
        case Experiment::Operation::BOOTSTRAP:
            return "BOOTSTRAP";
        case Experiment::Operation::ADD_SCALAR_SYNC:
            return "ADD_SCALAR_SYNC";
        case Experiment::Operation::ADD_SYNC:
            return "ADD_SYNC";
        case Experiment::Operation::ENCRYPT_SYNC:
            return "ENCRYPT_SYNC";
        case Experiment::Operation::MUL_SYNC:
            return "MUL_SYNC";
        case Experiment::Operation::MUL_SCALAR_SYNC:
            return "MUL_SCALAR_SYNC";
        case Experiment::Operation::ROT_SYNC:
            return "ROT_SYNC";
        case Experiment::Operation::ENCODE_SYNC:
            return "ENCODE_SYNC";
        case Experiment::Operation::CIPHER_SYNC:
            return "CIPHER_SYNC";
        case Experiment::Operation::ADD_SCALAR_COMPUTE:
            return "ADD_SCALAR_COMPUTE";
        case Experiment::Operation::ADD_COMPUTE:
            return "ADD_COMPUTE";
        case Experiment::Operation::MUL_COMPUTE:
            return "MUL_COMPUTE";
        case Experiment::Operation::MUL_SCALAR_COMPUTE:
            return "MUL_SCALAR_COMPUTE";
        case Experiment::Operation::ROT_COMPUTE:
            return "ROT_COMPUTE";
        case Experiment::Operation::ENCODE_COMPUTE:
            return "ENCODE_COMPUTE";
        case Experiment::Operation::ENCRYPT_COMPUTE:
            return "ENCRYPT_COMPUTE";
        case Experiment::Operation::ADD_SCALAR_PREP:
            return "ADD_SCALAR_PREP";
        case Experiment::Operation::ADD_PREP:
            return "ADD_PREP";
        case Experiment::Operation::MUL_PREP:
            return "MUL_PREP";
        case Experiment::Operation::MUL_SCALAR_PREP:
            return "MUL_SCALAR_PREP";
        case Experiment::Operation::ROT_PREP:
            return "ROT_PREP";
        case Experiment::Operation::ENCODE_PREP:
            return "ENCODE_PREP";
        case Experiment::Operation::ENCRYPT_PREP:
            return "ENCRYPT_PREP";
        case Experiment::Operation::PREP:
            return "PREP";
        case Experiment::Operation::MUL_PLAIN:
            return "MUL_PLAIN";
        case Experiment::Operation::MUL_PLAIN_SYNC:
            return "MUL_PLAIN_SYNC";
        case Experiment::Operation::MUL_PLAIN_COMPUTE:
            return "MUL_PLAIN_COMPUTE";
        case Experiment::Operation::MUL_PLAIN_PREP:
            return "MUL_PLAIN_PREP";

        default:
            return "UNKNOWN";
    }
}

std::string to_string(Experiment::Unit unit) {
    switch (unit) {
        case Experiment::Unit::MICROSECONDS:
            return "MICROSECONDS";
        case Experiment::Unit::MILLISECONDS:
            return "MILLISECONDS";
        case Experiment::Unit::SECONDS:
            return "SECONDS";
        case Experiment::Unit::UNIT:
            return "UNIT";
        default:
            return "UNKNOWN";
    }
}

Experiment::Experiment(std::string filename) {
    this->save_to_file = true;
    // Open the file
    this->file = new std::ofstream(filename, std::ios::out | std::ios::binary);
    if (!this->file->is_open()) {
#ifdef EXPE
        // Create folder expe
        mkdir(folder_expe.c_str(), 0777);
#endif
        this->file = new std::ofstream(filename, std::ios::out | std::ios::binary);
        if (!this->file->is_open()) {
            std::cerr << "Error: could not open file " << filename << std::endl;
            return;
        }
    }

    // Write the header
    *this->file << "multDepth,batchSize,ringsize,title,operation,unit,time,run" << std::endl;
}

void Experiment::saveToFile(std::string filename) {
    this->save_to_file = true;
    // Open the file
    this->file = new std::ofstream(filename, std::ios::out | std::ios::binary);
    if (!this->file->is_open()) {
#ifdef EXPE
        // Create folder expe
        mkdir(folder_expe.c_str(), 0777);
#endif
        this->file = new std::ofstream(filename, std::ios::out | std::ios::binary);
        if (!this->file->is_open()) {
            std::cerr << "Error: could not open file " << filename << std::endl;
            return;
        }
    }

    // Write the header
    *this->file << "multDepth,batchSize,ringsize,title,operation,unit,time,run" << std::endl;
}

void Experiment::save(std::string filename) {
    std::ofstream file(filename, std::ios::out | std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error: could not open file " << filename << std::endl;
        return;
    }

    // Serialize the current object
    file << "multDepth,batchSize,ringsize,title,operation,unit,time,run" << std::endl;
    for (auto statement : this->times) {
        if (crypto_context_available) {
            file << statement.crypto_params.multDepth << ",";
            file << statement.crypto_params.batchSize << ",";
            file << statement.ring_size << ",";
        } else {
            file << "N/A,N/A,N/A,";
        }
        file << statement.title << ",";
        file << to_string(statement.operation) << ",";
        file << to_string(statement.unit) << ",";
        file << statement.value << ",";
        file << statement.run << std::endl;
    }

    file.close();
}

void Experiment::setCryptoParams(CCParams<CryptoContextCKKSRNS> params) {
    this->current_params.multDepth = params.GetMultiplicativeDepth();
    this->current_params.batchSize = params.GetBatchSize();
    this->crypto_context_available = true;
}

void Experiment::setRingsize(int ringsize) {
    this->current_ring_size = ringsize;
    this->crypto_context_available = true;
}

void Experiment::addValue(std::vector<long> times, std::string title, Operation operation) {
    for (auto time : times) {
        this->times.push_back({this->current_params, time, title, operation, MICROSECONDS, this->current_ring_size});
        if (this->save_to_file) {
            if (crypto_context_available) {
                *this->file << this->current_params.multDepth << ",";
                *this->file << this->current_params.batchSize << ",";
                *this->file << this->current_ring_size << ",";
            } else {
                *this->file << "N/A,N/A,N/A,";
            }
            *this->file << title << ",";
            *this->file << to_string(operation) << ",";
            *this->file << to_string(MICROSECONDS) << ",";
            *this->file << time << ",";
            *this->file << this->current_run << std::endl;
        }
    }
}

void Experiment::addValue(long time, std::string title, Operation operation, Unit unit) {
    this->times.push_back({this->current_params, time, title, operation, unit, this->current_ring_size});
    if (this->save_to_file) {
        if (crypto_context_available) {
            *this->file << this->current_params.multDepth << ",";
            *this->file << this->current_params.batchSize << ",";
            *this->file << this->current_ring_size << ",";
        } else {
            *this->file << "N/A,N/A,N/A,";
        }
        *this->file << title << ",";
        *this->file << to_string(operation) << ",";
        *this->file << to_string(unit) << ",";
        *this->file << time << ",";
        *this->file << this->current_run << std::endl;
    }
}

void Experiment::addValue(std::vector<long> times, std::string title, std::vector<Operation> operations, Unit unit) {
    assert(times.size() == operations.size());
    for (int i = 0; i < times.size(); i++) {
        this->times.push_back({this->current_params, times[i], title, operations[i], unit, this->current_ring_size});
        if (this->save_to_file) {
            if (crypto_context_available) {
                *this->file << this->current_params.multDepth << ",";
                *this->file << this->current_params.batchSize << ",";
                *this->file << this->current_ring_size << ",";
            } else {
                *this->file << "N/A,N/A,N/A,";
            }
            *this->file << title << ",";
            *this->file << to_string(operations[i]) << ",";
            *this->file << to_string(unit) << ",";
            *this->file << times[i] << ",";
            *this->file << this->current_run << std::endl;
        }
    }
}

std::ostream& operator<<(std::ostream& os, const Experiment& dt) {
    os << "Experiment: " << std::endl;
    for (auto statement : dt.times) {
        os << "  CryptoParams: " << statement.crypto_params.multDepth << ", " << statement.crypto_params.batchSize
           << std::endl;
        os << "  Ringsize: " << statement.ring_size << std::endl;
        os << "  Title: " << statement.title << std::endl;
        os << "  Operation: " << statement.operation << std::endl;
        os << "  Unit: " << statement.unit << std::endl;
        os << "  Time: " << statement.value << std::endl;
    }
    return os;
}
const std::string currentDateTimeExpe() {
    time_t now = time(nullptr);
    struct tm tstruct;
    std::array<char, 80> buf;
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf.data(), sizeof(buf), "%Y-%m-%d.%X", &tstruct);

    return std::string(buf.data());
}