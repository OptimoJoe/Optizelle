// std::vector
#include <vector>

// std::unique_ptr
#include <memory>

// std::ctype
#include <locale>

namespace Optizelle { namespace Stream {
    // Basic stream type
    template <typename A>
    struct t {
        // Required functions
        virtual std::unique_ptr <A> next() = 0;

        // Allow for destructors
        virtual ~t() {}
    };

    // Converting a stream from the standard library
    template <typename Stream,typename A>
    struct of_std : public t <A> {
        // Internal stream
        std::unique_ptr <Stream> stream;

        // Mask for setting whitespace characters
        std::vector<std::ctype<char>::mask> mask;

        // Obtains memory for the stream and sets what delimiters we want for
        // whitespace.  Note, this is only really works for streams containing
        // characters, so use carefully.
        of_std(
            Stream * stream_,
            std::vector <char> delims = {'\n', ' ', '\t'}
        ) : stream(stream_), mask() {
            // Grab the normal table of characters
            auto const table = std::ctype<char>::classic_table();

            // Create a mask that we'll modify
            mask = std::vector<std::ctype<char>::mask>(
                table,table+std::ctype<char>::table_size);

            // Turn off the normal whitespace characters
            mask[' '] &= ~(std::ctype_base::space);
            mask['\t'] &= ~(std::ctype_base::space);
            mask['\n'] &= ~(std::ctype_base::space);

            // Turn on our new whitespace characters
            for(auto const & delim : delims)
                mask[delim] |= std::ctype_base::space;

            // Put this mask into the stream
            stream->imbue(
                std::locale(
                    stream->getloc(),new std::ctype<char>(mask.data())));
        }
        of_std(of_std &&) = default;
        of_std & operator = (of_std &&) = default;

        // Grab the next item in the stream
        std::unique_ptr <A> next() try {
            // Try to grab the next element
            auto item = std::unique_ptr <A> (new A());
            if(!(*stream >> *item))
                // Clear out what we return if we fail
                item.reset(nullptr);

            // Check for any really bad errors
            CHECK_STREAM(*stream);

            // Return the resulting element
            return item;
        } catch(...) {
            std::throw_with_nested(
                Exception::t(__LOC__
                    + ", problem grabbing the next element of the stream"));
        }
    };

    // Iterate over the stream
    template <typename A>
    void iter(std::function <void(A const &)> const & f,t <A> & stream) try {
        auto item = std::unique_ptr <A> (nullptr);
        while((item = stream.next()))
            f(*item);
    } catch(...) {
        std::throw_with_nested(
            Exception::t(__LOC__ +", problem iterating over the stream"));
    }

    // Iterate over the stream with a special item for the last element
    template <typename A>
    void iter_with_last(
        std::function <void(A const &)> const & f,
        std::function <void(A const &)> const & g,
        t <A> & stream
    ) try {
        // Attempt to read the next element and return if there is no element
        auto item = std::unique_ptr <A> (nullptr);
        if(!(item = stream.next()))
            return;

        // Since we're not empty, cache the element
        auto cache = std::move(item);

        // Continue looping while we're not empty
        while((item = stream.next())) {
            // Apply the function to the cached element
            f(*cache);

            // Move the next item into the cache
            cache = std::move(item);
        }

        // Once the stream is empty, then apply our second element to last
        // stream item
        g(*cache);
    } catch(...) {
        std::throw_with_nested(
            Exception::t(__LOC__ +", problem iterating over the stream "
                "with a special function for the last element"));
    }

    // Map a function to a stream
    template <typename A,typename B>
    struct map : public t <B> {
        // Function mapped to the stream
        std::function <B(A const &)> f;

        // Stream
        std::unique_ptr <typename Stream::t <A>> stream;

        // Grab the function mapped to the stream and the stream itself
        template <typename S>
        map(
            std::function <B(A const &)> const & f_,
            S && stream_
        ) : f(f_), stream(std::make_unique <S> (std::move(stream_))) {};

        // Grab the next item
        std::unique_ptr <B> next() {
            // If we can grab an element, apply the function to it and return
            auto item = std::unique_ptr <A> (nullptr);
            if((item = stream -> next()))
                return std::make_unique <B> (f(*item));

            // Otherwise, return null
            else
                return nullptr;
        }
    };

    // Filter elements from a stream
    template <typename A>
    struct filter : public t <A> {
        // Function to filter the stream
        std::function <bool(A const &)> const f;

        // Stream
        std::unique_ptr <t <A>> stream;

        // Grab the filter function and the stream itself
        template <typename S>
        filter(
            std::function <bool(A const &)> const & f_,
            S && stream_
        ) :
            f(f_),
            stream(std::make_unique <S> (std::move(stream_)))
        {}

        // Grab the next item
        std::unique_ptr <A> next() {
            // Loop over elements until we either find a valid one or we return
            // null
            auto item = std::unique_ptr <A> (nullptr);
            while((item = stream->next()) && !f(*item));

            // At this point, we either have a valid element or a null ptr.  In
            // either case, just return the result.
            return item;
        }
    };

    // Folds a stream into a function
    template <typename A,typename B>
    B fold(
        std::function<B(A const &,B const &)> const & f,
        t <A> & stream,
        B acc
    ) {
        // While we're not empty, keep grabbing elements
        auto item = std::unique_ptr <A> (nullptr);
        while((item = stream.next()))
            acc = f(*item,acc);

        // Return what we have
        return acc;
    }
}}
