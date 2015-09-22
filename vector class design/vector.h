// vector.h -- header file for vector data structure project


#include <cstring>
#include <stdexcept>

#pragma once
#ifndef _Vector_h
#define _Vector_h


namespace epl {
    
    template <typename T>
    class vector {
    private:
        T* buffer{nullptr};
        int64_t front{0};
        int64_t back{0};
        uint64_t capacity{0};
        uint64_t version{0};
        uint64_t assignment_version{0};
        uint64_t push_fronts{0};
        uint64_t pop_fronts{0};
        
    public:
        vector(void) {
            
            buffer = (T*) operator new(8 * sizeof(T));
            front = -1;
            back = -1;
            capacity = 8;
            
        }
        
        explicit vector(uint64_t n) {
            
            if (n == 0) {
                buffer = (T*) operator new(8 * sizeof(T));
                front = -1;
                back = -1;
                capacity = 8;
            } else {
                front = 0;
                back = n - 1;
                capacity = n;
                buffer = (T*) operator new(n * sizeof(T));
                for (int64_t i=0; i < n; i += 1) {
                    new (buffer + i) T{};
                }
            }
        }
        
        vector(const vector<T>& that) {
            
            copy(that);
            
        }
        
        vector(vector&& rhs) {
            
            this->front = rhs.front;
            this->back = rhs.back;
            this->capacity = rhs.capacity;
            this->buffer = rhs.buffer;
            
            rhs.buffer = nullptr;
            rhs.front = -1;
            rhs.back = -1;
            rhs.capacity = 0;
            rhs.version += 1;
            
        }
        
        vector<T>& operator=(const vector<T>& that) {
            if (this != &that) {
               
                for (int i=0; i < size(); ++i) {
                    (this->operator[](i)).~T();
                }
                
                if (that.capacity > this->capacity) {
                   
                    operator delete(buffer);
                    if (that.capacity == 0) {
                        this->buffer = nullptr;
                    } else {
                        this->buffer = (T*) operator new(that.capacity * sizeof(T));
                    }
                    this->capacity = that.capacity;
                }
                
                uint64_t that_size = that.size();
                
                if (that_size == 0) {
                    this->front = -1;
                    this->back = -1;
                } else {
                    this->front = 0;
                    this->back = that_size - 1;
                    
                    for(uint64_t i=0; i < that_size; i+=1) {
                        new (this->buffer + i) T{that[i]};
                    }
                }
                
                version += 1;
                assignment_version += 1;
            }
            return *this;
        }
        
        
        vector<T>& operator=(vector<T>&& rhs) {
            std::swap(this->buffer, rhs.buffer);
            std::swap(this->front, rhs.front);
            std::swap(this->back, rhs.back);
            std::swap(this->capacity, rhs.capacity);
            
            
            this->version += 1;
            rhs.version += 1;
            assignment_version += 1;
            
            return *this;
        }
        
        
        ~vector<T>(void) {
            destroy();
        }
        
        
        T& operator[](uint64_t index) {
            if (index >= size()) {
                throw std::out_of_range("subscript out of range");
            }
            return buffer[(front + index) % capacity];
        }
        
        
        const T& operator[](uint64_t index) const {
            if (index >= size()) {
                throw std::out_of_range("subscript out of range");
            }
            return buffer[(front + index) % capacity];
        }
        
        
        uint64_t size(void) const {
            if (front == -1) {
                return 0;
            }
            if (back >= front) {
                return back - front + 1;
            } else {
                return back + (capacity - front) + 1;
            }
        }
        void push_back(const T& arg) {
            
            if (front == (back + 1) % capacity) {
                
                T* new_loc = (T*) operator new(capacity * 2 * sizeof(T));
                
                
                new (new_loc + capacity) T(arg);
                
                for (int i=0; i < size(); ++i) {
                    new (new_loc + i) T{std::move(this->operator [](i))};
                }
                
                
                for (int i=0; i < size(); ++i) {
                    
                    (buffer + i)->~T();
                }
                operator delete(buffer);
                
                front = 0;
                back = capacity;
                capacity = capacity * 2;
                
                
                buffer = new_loc;
            } else {
                back = (back + 1) % capacity;
                new (buffer + back) T(arg);
                
                
                if (front == -1) {
                    front = 0;
                }
            }
            
            
            version += 1;
        }
        
        
        void push_back(T&& arg) {
            
            if (front == (back + 1) % capacity) {
                
                T* new_loc = (T*) operator new(capacity * 2 * sizeof(T));
                
                for (int i=0; i < size(); ++i) {
                    new (new_loc + i) T{std::move(this->operator [](i))};
                }
                
                new (new_loc + capacity) T(std::move(arg));
                
                
                for (int i=0; i < size(); ++i) {
                    
                    (buffer + i)->~T();
                }
                operator delete(buffer);
                
                front = 0;
                back = capacity;
                capacity = capacity * 2;
                
                buffer = new_loc;
            } else {
                back = (back + 1) % capacity;
                
                new (buffer + back) T(std::move(arg));
                
                
                if (front == -1) {
                    front = 0;
                }
            }
            
            version += 1;
        }
        
        
        void push_front(const T& arg) {
            
            if (front == (back + 1) % capacity) {
              
                T* new_loc = (T*) operator new(capacity * 2 * sizeof(T));
                
                
                new (new_loc + (capacity * 2) - 1) T(arg);
                
                for (int i=0; i < size(); ++i) {
                    new (new_loc + i) T{std::move(this->operator [](i))};
                }
                
                
                for (int i=0; i < size(); ++i) {
                    (buffer + i)->~T();
                }
                operator delete(buffer);
                
                back = capacity - 1;
                capacity = capacity * 2;
                front = capacity - 1;
                
                buffer = new_loc;
            } else {
                if (front == -1) {
                    front = 0;
                    back = 0;
                } else {
                    if (front == 0) {
                        front = capacity;
                    }
                    front = front - 1;
                }
                new (buffer + front) T(arg);
            }
            
            version += 1;
            push_fronts += 1;
        }
        
        void push_front(T&& arg) {
            
            if (front == (back + 1) % capacity) {
                
                T* new_loc = (T*) operator new(capacity * 2 * sizeof(T));
                
                for (int i=0; i < size(); ++i) {
                    new (new_loc + i) T{std::move(this->operator [](i))};
                }
                
                
                new (new_loc + (capacity * 2) - 1) T(std::move(arg));
                
                for (int i=0; i < size(); ++i) {
                    
                    (buffer + i)->~T();
                }
                operator delete(buffer);
                
                back = capacity - 1;
                capacity = capacity * 2;
                front = capacity - 1;
                buffer = new_loc;
            } else {
                if (front == -1) {
                    front = 0;
                    back = 0;
                } else {
                    if (front == 0) {
                        front = capacity;
                    }
                    front = front - 1;
                }
                
                new (buffer + front) T(std::move(arg));
            }
            
            version += 1;
            push_fronts += 1;
        }
        
        void pop_back(void) {
            
            if (front == -1) {
                throw std::out_of_range("vector is already empty.");
            }
            
            (buffer + back)->~T();
            
            if (back == 0) {
                back = capacity;
            }
            back -= 1;
            
            if (front == (back + 1) % capacity) {
                back = -1;
                front = -1;
            }

            version += 1;
        }
        
        void pop_front(void) {
            if (front == -1) {
                throw std::out_of_range("vector is already empty.");
            }
            
            (buffer + front)->~T();
            
            front = (front + 1) % capacity;

            if (front == (back + 1) % capacity) {
                back = -1;
                front = -1;
            }
            
            version += 1;
            pop_fronts += 1;
            
        }
        
    private:
        
        void destroy(void) {
            for (int i=0; i < size(); ++i) {
                (this->operator [](i)).~T();
            }
            operator delete(buffer);
        }
        
        void copy(const vector<T>& that) {
            this->front = that.front;
            this->back = that.back;
            this->capacity = that.capacity;
            
            if (capacity == 0) {
                this->buffer = nullptr;
            }
            else
            {
                this->buffer = (T*) operator new(capacity * sizeof(T));
                for (uint64_t k=0; k < size(); k += 1) {
                    new (this->buffer + ((front + k) % capacity)) T{that[k]};
                }
            }
        }
    };
    
} //namespace epl

#endif /* _vector_h */