#pragma once

#include <functional>
#include <type_traits>
#include <utility>
#include <variant>

namespace result_detail {
template <class... Ts>
struct overloaded : Ts... {
    using Ts::operator()...;
};

template <class... Ts>
overloaded(Ts...) -> overloaded<Ts...>;
} // namespace result_detail

template <class T>
struct AutoCastFailType {
    std::decay_t<T> value;
};

template <class T>
inline auto return_failure(T&& t) {
    return AutoCastFailType<T> { std::forward<T>(t) };
}

/*!
 * \brief Small success-or-failure value type.
 *
 * Result stores either a Success value or a Failure value. It is intended for
 * operations that can fail with a useful error object, while keeping success
 * values explicit at call sites.
 *
 * This type inherits std::variant so existing variant construction remains
 * available. Success and Failure must be different types because the active
 * state is selected by type.
 */
template <class Success, class Failure>
class Result : public std::variant<Success, Failure> {
    static_assert(!std::is_same_v<Success, Failure>,
                  "Same class for success and failure are not supported");

public:
    using Base = std::variant<Success, Failure>;

    /*!
     * \brief Default-constructs the underlying variant.
     *
     * Like std::variant<Success, Failure>, this constructs the first
     * alternative, so the default state is Success.
     */
    Result() = default;

    template <class T>
    Result(AutoCastFailType<T> ftype)
        : Result(Failure(std::move(ftype.value))) { }

    /// Reuse std::variant constructors for Success and Failure values.
    using Base::Base;

    /// Returns true when this object currently contains a Success.
    bool is_success() const {
        return std::get_if<Success>(base_ptr()) != nullptr;
    }

    /// Returns true when this object currently contains a Failure.
    bool is_failure() const {
        return std::get_if<Failure>(base_ptr()) != nullptr;
    }

    /// Returns the contained Success, or throws std::bad_variant_access.
    Success& get_success() { return std::get<Success>(base_ref()); }

    /// Returns the contained Failure, or throws std::bad_variant_access.
    Failure& get_failure() { return std::get<Failure>(base_ref()); }

    /// Returns the contained Success, or throws std::bad_variant_access.
    Success const& get_success() const { return std::get<Success>(base_ref()); }

    /// Returns the contained Failure, or throws std::bad_variant_access.
    Failure const& get_failure() const { return std::get<Failure>(base_ref()); }

    /*!
     * \brief Calls one of two functions with a const reference to the value.
     *
     * The success function is called for Success, and the failure function is
     * called for Failure. The contained value is not moved.
     */
    template <class SFunction, class FFunction>
    decltype(auto) match_ref(SFunction&& sf, FFunction&& ff) const {
        return std::visit(
            result_detail::overloaded {
                [&sf](Success const& s) -> decltype(auto) {
                    return std::invoke(std::forward<SFunction>(sf), s);
                },
                [&ff](Failure const& f) -> decltype(auto) {
                    return std::invoke(std::forward<FFunction>(ff), f);
                },
            },
            static_cast<Base const&>(*this));
    }

    /*!
     * \brief Calls one of two functions while consuming the contained value.
     *
     * This overload is rvalue-qualified. Use std::move(result).match_consume(...)
     * when the handler should receive Success&& or Failure&&.
     */
    template <class SFunction, class FFunction>
    decltype(auto) match_consume(SFunction&& sf, FFunction&& ff) && {
        return std::visit(
            result_detail::overloaded {
                [&sf](Success&& s) -> decltype(auto) {
                    return std::invoke(std::forward<SFunction>(sf),
                                       std::move(s));
                },
                [&ff](Failure&& f) -> decltype(auto) {
                    return std::invoke(std::forward<FFunction>(ff),
                                       std::move(f));
                },
            },
            static_cast<Base&&>(std::move(*this)));
    }

    /*!
     * \brief Transforms a Success value while preserving Failure.
     *
     * For lvalue Results, Success is passed to the transform as a const
     * reference and Failure is copied into the returned Result.
     */
    template <class NewSuccess, class Function>
    Result<NewSuccess, Failure> map(Function&& f) const& {
        return match_ref(
            [&f](Success const& x) {
                return Result<NewSuccess, Failure>(
                    std::invoke(std::forward<Function>(f), x));
            },
            [](Failure const& x) { return Result<NewSuccess, Failure>(x); });
    }

    /*!
     * \brief Transforms a Success value while consuming this Result.
     *
     * For rvalue Results, Success and Failure are moved into the transform or
     * returned Result.
     */
    template <class NewSuccess, class Function>
    Result<NewSuccess, Failure> map(Function&& f) && {
        return std::move(*this).match_consume(
            [&f](Success&& x) {
                return Result<NewSuccess, Failure>(
                    std::invoke(std::forward<Function>(f), std::move(x)));
            },
            [](Failure&& x) {
                return Result<NewSuccess, Failure>(std::move(x));
            });
    }

    /*!
     * \brief Transforms a Failure value while preserving Success.
     *
     * For lvalue Results, Failure is passed to the transform as a const
     * reference and Success is copied into the returned Result.
     */
    template <class NewFailure, class Function>
    Result<Success, NewFailure> map_error(Function&& f) const& {
        return match_ref(
            [](Success const& x) { return Result<Success, NewFailure>(x); },
            [&f](Failure const& x) {
                return Result<Success, NewFailure>(
                    std::invoke(std::forward<Function>(f), x));
            });
    }

    /*!
     * \brief Transforms a Failure value while consuming this Result.
     *
     * For rvalue Results, Success and Failure are moved into the transform or
     * returned Result.
     */
    template <class NewFailure, class Function>
    Result<Success, NewFailure> map_error(Function&& f) && {
        return std::move(*this).match_consume(
            [](Success&& x) {
                return Result<Success, NewFailure>(std::move(x));
            },
            [&f](Failure&& x) {
                return Result<Success, NewFailure>(
                    std::invoke(std::forward<Function>(f), std::move(x)));
            });
    }

    /*!
     * \brief Converts to another Result when both alternatives are convertible.
     */
    template <class OtherSuccess, class OtherError>
    operator Result<OtherSuccess, OtherError>() const {
        static_assert(std::is_convertible_v<Success, OtherSuccess>);
        static_assert(std::is_convertible_v<Failure, OtherError>);

        return map([](auto x) { return OtherSuccess(x); })
            .map_error([](auto x) { return OtherError(x); });
    }

    /// Shorthand for is_success().
    operator bool() const { return is_success(); }

private:
    Base&       base_ref() { return static_cast<Base&>(*this); }
    Base const& base_ref() const { return static_cast<Base const&>(*this); }
    Base const* base_ptr() const { return &base_ref(); }
};
