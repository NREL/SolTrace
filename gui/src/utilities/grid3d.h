#pragma once

#include <array>
#include <span>
#include <stdexcept>

#include <QHash>
#include <QVector>

#define GLM_ENABLE_EXPERIMENTAL 1
#include <glm/common.hpp>
#include <glm/gtx/hash.hpp>
#include <glm/vec3.hpp>

namespace analysis {

/*!
 * \brief A 3D grid, with integral indicies.
 */

template <class T>
class Grid3D {
    QVector<T> m_data;

    using DimType = std::array<size_t, 3>;
    DimType m_dimensions;

    inline void init_size() {
        m_data.resize(m_dimensions[0] * m_dimensions[1] * m_dimensions[2]);
    }

public:
    /*!
     * \brief Construct an empty #D grid (x and y dimensions are 0).
     */
    explicit Grid3D() : Grid3D(0, 0, 0) { }

    /*!
     * \brief Construct a 3d grid with the specified dimensions.
     */
    explicit Grid3D(DimType const& dims) : m_dimensions(dims) { init_size(); }

    /*!
     * \brief Construct a 3d grid with the specified dimensions.
     */
    explicit Grid3D(size_t xd, size_t yd, size_t zd)
        : m_dimensions({ { xd, yd, zd } }) {
        init_size();
    }

    ~Grid3D() = default;

    /*!
     * \brief Set all grid points to the given value.
     */
    inline void fill(T const& t) {
        for (T& lt : m_data)
            lt = t;
    }

    /*!
     * \brief The number of grid points.
     */
    inline unsigned size() const { return m_data.size(); }

    /*!
     * \brief Compute the linear offset of a given point.
     *
     * This can be used to index into an array
     */
    inline unsigned index(size_t x, size_t y, size_t z) const {
        return x + m_dimensions[0] * (y + m_dimensions[1] * z);
    }

    /*!
     * \brief The x dimension.
     */
    inline size_t size_x() const { return m_dimensions[0]; }

    /*!
     * \brief The y dimension.
     */
    inline size_t size_y() const { return m_dimensions[1]; }

    /*!
     * \brief The z dimension.
     */
    inline size_t size_z() const { return m_dimensions[2]; }


    /*!
     * \brief Access a grid point.
     */
    T& operator()(size_t x, size_t y, size_t z) {
        return m_data[index(x, y, z)];
    }

    /*!
     * \brief Access a grid point.
     */
    T const& operator()(size_t x, size_t y, size_t z) const {
        return m_data[index(x, y, z)];
    }

    /*!
     * \brief Access a grid point by its linear index.
     *
     * Useful for iterating through the grid.
     */
    T& operator[](size_t x) { return m_data[x]; }

    /*!
     * \brief Access a grid point by its linear index.
     *
     * Useful for iterating through the grid.
     */
    T const& operator[](size_t x) const { return m_data[x]; }

    /*!
     * \brief Access a grid point.
     *
     * This will throw if out of bounds.
     */
    T& at(size_t x, size_t y, size_t z) {
        auto i = index(x, y, z);
        if (i >= size()) throw std::out_of_range("Grid index out of range");
        return m_data[i];
    }

    /*!
     * \brief Access a grid point.
     *
     * This will throw if out of bounds.
     */
    T const& at(size_t x, size_t y, size_t z) const {
        auto i = index(x, y, z);
        if (i >= size()) throw std::out_of_range("Grid index out of range");
        return m_data[i];
    }

    /*!
     * \brief Ensure that two given integers are within the dimensions of the
     * grid.
     *
     * Useful to make sure accesses are not out of bound.
     */
    void clamp_bounds(size_t& x, size_t& y, size_t& z) const {
        x = clamp(x, 0u, size_x() - 1);
        y = clamp(y, 0u, size_y() - 1);
        z = clamp(z, 0u, size_z() - 1);
    }

public:
    /*!
     * \brief Linear iterator support.
     */
    auto begin() { return m_data.begin(); }

    /*!
     * \brief Linear iterator support.
     */
    auto end() { return m_data.end(); }

    /*!
     * \brief Access the underlying storage.
     */
    std::span<T>& as_vector() { return m_data; }

    /*!
     * \brief Access the underlying storage.
     */
    std::span<T> const& as_vector() const { return m_data; }

    /*!
     * \brief Get a pointer to the first element in the grid.
     */
    T*       data() { return m_data.data(); }
    T const* data() const { return m_data.data(); }
};


/*!
 * \brief A 3D grid, with integral indicies.
 */

template <class T>
class SparseGrid3D {
public:
    template <class U>
    struct Brick {
        std::array<U, 16 * 16 * 16> values = {};

        T& operator()(glm::ivec3 p) {
            auto at = p.x + 16 * (p.y + 16 * p.z);
            return values[at];
        }

        T const& operator()(glm::ivec3 p) const {
            auto at = p.x + 16 * (p.y + 16 * p.z);
            return values[at];
        }
    };


    explicit SparseGrid3D(T background = T()) : m_background(background) { }

    SparseGrid3D(SparseGrid3D&&)            = default;
    SparseGrid3D& operator=(SparseGrid3D&&) = default;

    SparseGrid3D(SparseGrid3D const&)            = default;
    SparseGrid3D& operator=(SparseGrid3D const&) = default;

    ~SparseGrid3D() = default;

    void update_active_bounds() {
        if (m_data.size() == 0) {
            m_active_min = glm::ivec3 {};
            m_active_max = glm::ivec3 {};
            return;
        }

        m_active_min = m_active_max = m_data.asKeyValueRange().begin()->first;

        for (auto const& [k, v] : m_data.asKeyValueRange()) {
            m_active_min = glm::min(m_active_min, k);
            m_active_max = glm::max(m_active_max, k);
        }

        // these are brick coords. need item coords?
        m_active_min *= 16;
        m_active_max = (m_active_max + 1) * 16;
    }

    T const& background() const { return m_background; }

    void set_transform(glm::vec3 translate, glm::vec3 scale) {
        m_translate = translate;
        m_scale     = scale;
    }

    glm::vec3 translate() const { return m_translate; }

    glm::vec3 scale() const { return m_scale; }

    glm::vec3 grid_to_world(glm::vec3 grid_position) const {
        return (grid_position - m_translate) / m_scale;
    }

    glm::vec3 world_to_grid(glm::vec3 world_position) const {
        return world_position * m_scale + m_translate;
    }

    glm::ivec3 active_min() const { return m_active_min; }
    glm::ivec3 active_max() const { return m_active_max; }

    glm::ivec3 active_span() const { return (m_active_max - m_active_min); }

    static inline int floor_div_16(int value) {
        if (value >= 0) { return value / 16; }

        return -(((-value) + 15) / 16);
    }

    static inline std::pair<glm::ivec3, glm::ivec3>
    block_brick_coord(glm::ivec3 p) {
        auto brick = glm::ivec3 {
            floor_div_16(p.x),
            floor_div_16(p.y),
            floor_div_16(p.z),
        };

        auto block = glm::ivec3 {
            p.x - brick.x * 16,
            p.y - brick.y * 16,
            p.z - brick.z * 16,
        };

        return { brick, block };
    }

    /*!
     * \brief The number of grid points.
     */
    inline unsigned size_in_bricks() const { return m_data.size(); }

    /*!
     * \brief Access a grid point.
     */
    T& operator()(glm::ivec3 p) {
        auto [brick, block] = block_brick_coord(p);
        return m_data[brick](block);
    }

    /*!
     * \brief Access a grid point.
     */
    T& operator()(int x, int y, int z) { return (*this)({ x, y, z }); }

    /*!
     * \brief Access a grid point.
     */
    T const& operator()(glm::ivec3 p) const {

        auto [brick, block] = block_brick_coord(p);

        auto iter = m_data.find(brick);

        if (iter == m_data.end()) { return m_background; }

        return iter.value()(block);
    }

    T const& operator()(int x, int y, int z) const {
        return (*this)({ x, y, z });
    }

    QHash<glm::ivec3, Brick<T>> const& data() const { return m_data; }

public:
    /*!
     * \brief Linear iterator support.
     */
    auto begin() { return m_data.begin(); }

    /*!
     * \brief Linear iterator support.
     */
    auto end() { return m_data.end(); }

private:
    // not the best...very slow of course
    T                           m_background = T();
    QHash<glm::ivec3, Brick<T>> m_data;

    // for perf reasons, we do not update these online
    glm::ivec3 m_active_min = {};
    glm::ivec3 m_active_max = {};

    glm::vec3 m_translate = glm::vec3(0.0f);
    glm::vec3 m_scale     = glm::vec3(1.0f);
};

template <class T>
class SparseGrid3DLookupCache {
    using GT    = SparseGrid3D<T>;
    using Brick = GT::template Brick<T>;

    GT const&                       m_source_grid;
    QHash<glm::ivec3, Brick> const& m_source;

    mutable glm::ivec3   m_current_brick_id = glm::ivec3(INT_MAX);
    mutable Brick const* m_current_brick    = nullptr;

public:
    SparseGrid3DLookupCache(SparseGrid3D<T> const& grid)
        : m_source_grid(grid), m_source(grid.data()) { }

    T const& operator()(glm::ivec3 p) const {
        auto [brick, block] = GT::block_brick_coord(p);

        if (m_current_brick_id != brick) {
            m_current_brick_id = brick;

            auto iter = m_source.find(brick);

            if (iter == m_source.end()) {
                m_current_brick = nullptr;
            } else {
                m_current_brick = &(iter.value());
            }
        }

        if (m_current_brick) { return (*m_current_brick)(block); }

        return m_source_grid.background();
    }

    T const& operator()(int x, int y, int z) const {
        return (*this)({ x, y, z });
    }
};

} // namespace analysis
