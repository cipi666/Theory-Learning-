# 算子学习的一些数学基础

## 奇异值分解

### 一、 SVD 的核心定义

对于任意一个 $m \times n$ 的实数矩阵 $A$，SVD 将其分解为：

$$A = U \Sigma V^T$$

这三个矩阵的含义如下：

1.  **$U$ (左奇异向量矩阵)**：
    * 大小：$m \times m$
    * 性质：是一个**正交矩阵** ($U^T U = I$)。它的列向量被称为**左奇异向量**。
    * 几何意义：代表在输出空间（Range）的旋转或反射。
2.  **$\Sigma$ (奇异值矩阵)**：
    * 大小：$m \times n$ (与 $A$ 形状相同，但在计算中常只写出非零部分)
    * 性质：是一个**对角矩阵**（除对角线外全为0）。对角线上的元素 $\sigma_i$ 称为**奇异值**。
    * 排序：通常按从大到小排列 ($\sigma_1 \geq \sigma_2 \geq \dots \geq 0$)。
    * 几何意义：代表沿坐标轴的**缩放**（拉伸或压缩）。
3.  **$V^T$ (右奇异向量矩阵的转置)**：
    * 大小：$n \times n$
    * 性质：$V$ 是**正交矩阵**。它的列向量被称为**右奇异向量**。
    * 几何意义：代表在输入空间（Domain）的旋转或反射。

**一句话总结几何意义**：任何线性变换 $A$，都可以看作是：先旋转 ($V^T$)，再沿坐标轴缩放 ($\Sigma$)，最后再旋转 ($U$)。

### 二、 计算步骤 (Step-by-Step)

为了手动计算 SVD，我们需要利用 $A$ 与其转置矩阵 $A^T$ 的关系。

#### 核心原理：
* $A^T A$ 是一个对称实矩阵，其特征向量组成了 $V$。
* $A A^T$ 也是一个对称实矩阵，其特征向量组成了 $U$。
* $A^T A$ 和 $A A^T$ 拥有相同的**非零特征值**，这些特征值的平方根就是**奇异值** $\sigma$。

#### 具体步骤：

假设我们有一个矩阵 $A$。

**第一步：求右奇异向量 $V$ 和奇异值 $\sigma$**
1.  计算方阵 $W = A^T A$。
2.  求 $W$ 的**特征值** $\lambda_i$ 和对应的**特征向量** $v_i$。
3.  将特征值从大到小排序：$\lambda_1 \geq \lambda_2 \dots$。
4.  **计算奇异值**：$\sigma_i = \sqrt{\lambda_i}$。
5.  **归一化**特征向量 $v_i$，得到单位向量，这些就是 $V$ 的列向量。

**第二步：求左奇异向量 $U$**
虽然我们可以通过计算 $A A^T$ 的特征向量来求 $U$，但这通常更麻烦且容易出现正负号不匹配的问题。**更简便且标准的方法**利用了 SVD 的定义 $Av_i = \sigma_i u_i$：

$$u_i = \frac{1}{\sigma_i} A v_i$$

1.  对于每一个非零的奇异值 $\sigma_i$，直接用上述公式计算 $u_i$。
2.  **注意**：如果 $m > n$（矩阵是竖长的），或者有零奇异值，上述公式求出的 $u$ 数量不够填满 $m \times m$ 的矩阵。此时，需要利用**施密特正交化 (Gram-Schmidt)** 找到剩余的正交向量来补全 $U$。

**第三步：组装**
将算出的 $U, \Sigma, V^T$ 组合起来。

### 三、 手算实例

让我们通过一个简单的例子来演示。
设矩阵 $A = \begin{pmatrix} 1 & 1 \\ 2 & 2 \end{pmatrix}$。这是一个 $2 \times 2$ 矩阵，且显然秩为 1（第二行是第一行的2倍）。

#### 1. 计算 $A^T A$
$$A^T A = \begin{pmatrix} 1 & 2 \\ 1 & 2 \end{pmatrix} \begin{pmatrix} 1 & 1 \\ 2 & 2 \end{pmatrix} = \begin{pmatrix} 1+4 & 1+4 \\ 1+4 & 1+4 \end{pmatrix} = \begin{pmatrix} 5 & 5 \\ 5 & 5 \end{pmatrix}$$

#### 2. 求 $A^T A$ 的特征值和特征向量
求解特征方程 $|A^T A - \lambda I| = 0$:
$$\begin{vmatrix} 5-\lambda & 5 \\ 5 & 5-\lambda \end{vmatrix} = (5-\lambda)^2 - 25 = \lambda^2 - 10\lambda = 0$$
$$\lambda(\lambda - 10) = 0$$
得到特征值：$\lambda_1 = 10, \quad \lambda_2 = 0$。

**奇异值 $\Sigma$**：
$\sigma_1 = \sqrt{10}, \quad \sigma_2 = \sqrt{0} = 0$。
所以 $\Sigma = \begin{pmatrix} \sqrt{10} & 0 \\ 0 & 0 \end{pmatrix}$。

**对应特征向量 ($V$)**：
* 对于 $\lambda_1 = 10$:
    $\begin{pmatrix} -5 & 5 \\ 5 & -5 \end{pmatrix} \begin{pmatrix} x \\ y \end{pmatrix} = 0 \Rightarrow x = y$。
    取单位向量：$v_1 = \begin{pmatrix} \frac{1}{\sqrt{2}} \\ \frac{1}{\sqrt{2}} \end{pmatrix}$。
* 对于 $\lambda_2 = 0$:
    $\begin{pmatrix} 5 & 5 \\ 5 & 5 \end{pmatrix} \begin{pmatrix} x \\ y \end{pmatrix} = 0 \Rightarrow x = -y$。
    取单位向量：$v_2 = \begin{pmatrix} \frac{1}{\sqrt{2}} \\ -\frac{1}{\sqrt{2}} \end{pmatrix}$。

所以 $V = \begin{pmatrix} \frac{1}{\sqrt{2}} & \frac{1}{\sqrt{2}} \\ \frac{1}{\sqrt{2}} & -\frac{1}{\sqrt{2}} \end{pmatrix}$。

#### 3. 求左奇异向量 $U$
利用公式 $u_i = \frac{1}{\sigma_i} A v_i$。

* 对于 $i=1$ ($\sigma_1 = \sqrt{10}$):
    $$u_1 = \frac{1}{\sqrt{10}} \begin{pmatrix} 1 & 1 \\ 2 & 2 \end{pmatrix} \begin{pmatrix} \frac{1}{\sqrt{2}} \\ \frac{1}{\sqrt{2}} \end{pmatrix} = \frac{1}{\sqrt{10}} \begin{pmatrix} \frac{2}{\sqrt{2}} \\ \frac{4}{\sqrt{2}} \end{pmatrix} = \begin{pmatrix} \frac{2}{\sqrt{20}} \\ \frac{4}{\sqrt{20}} \end{pmatrix} = \begin{pmatrix} \frac{1}{\sqrt{5}} \\ \frac{2}{\sqrt{5}} \end{pmatrix}$$

* 对于 $i=2$ ($\sigma_2 = 0$):
    无法使用公式（分母为0）。我们需要找一个与 $u_1$ 正交的单位向量。
    $u_1 = \begin{pmatrix} 1/\sqrt{5} \\ 2/\sqrt{5} \end{pmatrix}$，与其正交的向量可以是 $\begin{pmatrix} -2/\sqrt{5} \\ 1/\sqrt{5} \end{pmatrix}$ 或 $\begin{pmatrix} 2/\sqrt{5} \\ -1/\sqrt{5} \end{pmatrix}$。
    这里我们要构建标准正交基，任选一个即可，设 $u_2 = \begin{pmatrix} -2/\sqrt{5} \\ 1/\sqrt{5} \end{pmatrix}$。

所以 $U = \begin{pmatrix} \frac{1}{\sqrt{5}} & -\frac{2}{\sqrt{5}} \\ \frac{2}{\sqrt{5}} & \frac{1}{\sqrt{5}} \end{pmatrix}$。

#### 4. 最终结果
$$A = U \Sigma V^T = \begin{pmatrix} \frac{1}{\sqrt{5}} & -\frac{2}{\sqrt{5}} \\ \frac{2}{\sqrt{5}} & \frac{1}{\sqrt{5}} \end{pmatrix} \begin{pmatrix} \sqrt{10} & 0 \\ 0 & 0 \end{pmatrix} \begin{pmatrix} \frac{1}{\sqrt{2}} & \frac{1}{\sqrt{2}} \\ \frac{1}{\sqrt{2}} & -\frac{1}{\sqrt{2}} \end{pmatrix}$$
*(注：$V^T$ 是 $V$ 的转置)*

### 四、 SVD 在视觉任务中的直观作用

既然你之前问过图像重建，SVD 在这里有一个非常关键的应用：**低秩近似 (Low-Rank Approximation)**。

如果我们将 $\Sigma$ 中较小的奇异值强行置为 0，只保留最大的 $k$ 个奇异值，再乘回去：
$$A_{approx} = U_k \Sigma_k V_k^T$$

这就实现了：
1.  **图像压缩**：用更少的数据存储图像。
2.  **去噪**：噪声通常对应于那些微小的奇异值，扔掉它们等于扔掉噪声。
3.  **特征提取**：最大的几个奇异值代表了图像的主要结构（比如人脸的轮廓），这就是主成分分析 (PCA) 的数学基础。


## 循环矩阵与傅里叶视角
在深度学习（特别是 FNO）和信号处理中，我们经常遇到巨大的矩阵运算。如果直接计算，成本太高。
本讲义的核心目标是解释一句话：
> **“在时域（空间域）中复杂的循环卷积，在傅里叶域中只是简单的逐元素乘法。”**

### 什么是循环矩阵 (Circulant Matrix)？

#### 定义
循环矩阵是一种非常特殊的 Toeplitz 矩阵。它的每一行都是上一行向右循环移动一个位置得到的。
由向量 $c = [c_0, c_1, \dots, c_{N-1}]^T$ 生成的循环矩阵 $C_c$ 如下所示：

$$
C_c = \begin{bmatrix} 
c_0 & c_{N-1} & \cdots & c_2 & c_1 \\
c_1 & c_0 & c_{N-1} & \ddots & c_2 \\
\vdots & c_1 & c_0 & \ddots & \vdots \\
c_{N-2} & \ddots & \ddots & \ddots & c_{N-1} \\
c_{N-1} & c_{N-2} & \cdots & c_1 & c_0 
\end{bmatrix}
$$

#### 物理意义：卷积
循环矩阵最重要的性质是：**矩阵与向量的乘法等价于卷积运算。**
假设有一个输入信号向量 $x$，计算 $y = C_c x$。
这在数学上完全等价于向量 $c$ 和 $x$ 的**循环卷积 (Circular Convolution)**：
$$y = c * x$$

* **直观理解**：矩阵的每一行就像是卷积核 $c$ 在信号 $x$ 上滑动的一个“快照”。

---

### 傅里叶坐标系

要理解“傅里叶基下的乘法算子”，我们首先要改变看待向量的方式。

#### 换个坐标系看世界
* **标准基（Standard Basis / Pixel Domain）：**
    当我们写一个向量 $x = [1, 2, 3]$ 时，我们默认使用的是标准基。意思是：第1个位置强度是1，第2个是2... 这是我们在**空间域**看到的图像。
    
* **傅里叶基（Fourier Basis / Frequency Domain）：**
    傅里叶变换告诉我们：任何信号都可以分解为不同频率的正弦/余弦波的叠加。
    **傅里叶基**就是一组特定频率的波。
    
    当我们做离散傅里叶变换 (DFT) 时，我们实际上是在做**基底变换 (Change of Basis)**。我们将向量从“像素坐标系”旋转到了“频率坐标系”。

#### 离散傅里叶变换矩阵 ($F$)
这种坐标变换可以通过一个矩阵 $F$（DFT矩阵）来实现：
* **正变换**（进入频率坐标系）：$\hat{x} = F x$
* **逆变换**（回到标准坐标系）：$x = F^{-1} \hat{x}$

### 核心定理：卷积定理与对角化

为什么我们说 $C_c$ 在傅里叶基下是“乘法算子”？

#### 卷积定理
数学上有一个著名的定理：
> **时域的卷积 等于 频域的逐元素相乘 (Pointwise Multiplication)。**

公式表达为：
$$\text{DFT}(c * x) = \text{DFT}(c) \odot \text{DFT}(x)$$
其中 $\odot$ 代表哈达玛积（Hadamard product），即对应位置元素直接相乘。

#### 矩阵语言的翻译（对角化）
如果我们把上面的定理翻译成线性代数的矩阵语言：
1.  **在原坐标系下**：运算是 $y = C_c x$ （这是稠密矩阵乘法）。
2.  **换到傅里叶坐标系**：
    * 先变身：把 $x$ 变成 $\hat{x} = Fx$。
    * 做运算：卷积变成了“逐元素乘法”。**在线性代数中，逐元素乘法向量 $\hat{c}$ 等价于乘以对角矩阵 $\text{diag}(\hat{c})$。**
    * 变回来：$\hat{y}$ 变回 $y = F^{-1} \hat{y}$。

把这三步连起来，我们得到：
$$C_c = F^{-1} \cdot \text{diag}(\hat{c}) \cdot F$$

这说明了：**循环矩阵 $C_c$ 可以被傅里叶矩阵 $F$ 对角化！**

* **特征向量**：$F$ 的列（即傅里叶基向量）就是 $C_c$ 的特征向量。
* **特征值**：$\hat{c}$ （即向量 $c$ 的傅里叶变换）就是对应的特征值。

#### 什么是“傅里叶基下的乘法算子”？
现在回顾那句难懂的话：
> "$C_c$ can be interpreted as a multiplication operator in the Fourier basis."

它的意思是：如果你站在傅里叶坐标系里看，$C_c$ 不再是一个乱糟糟的满矩阵，它变成了一个干净的**对角矩阵**。对角矩阵的作用就是把向量的每个分量单独乘以一个数（缩放），这就是所谓的“乘法算子”。

### 为什么这很重要？（计算复杂度）

理解了这个结构，我们就能明白为什么 FNO 和相关算法如此高效。

* **直接计算 $y = C_c x$**：
    $C_c$ 是 $N \times N$ 矩阵。
    计算复杂度：**$O(N^2)$**。
    
* **利用傅里叶技巧计算**：
    利用分解 $C_c x = F^{-1} (\text{diag}(\hat{c}) (F x))$。
    1.  FFT ($F x$)：$O(N \log N)$
    2.  逐元素乘法 ($\odot \hat{c}$)：$O(N)$
    3.  IFFT ($F^{-1}$)：$O(N \log N)$
    
    总复杂度：**$O(N \log N)$**。
    
    **结论**：当 $N$ 很大时（例如图像像素），$N \log N$ 比 $N^2$ 快几个数量级。

### 附录：关于 $C_c g = C_g c$ 的解释

您提供的截图中提到了利用 $g$ 恢复 $c$ 的技巧，利用了卷积的**交换律**。

* **卷积视角**：$c * g = g * c$
* **矩阵视角**：
    * $c * g$ 对应矩阵运算 $C_c g$ （固定 $c$，把 $c$ 做成矩阵去卷 $g$）。
    * $g * c$ 对应矩阵运算 $C_g c$ （固定 $g$，把 $g$ 做成矩阵去卷 $c$）。
    
因为卷积结果一样，所以这两个矩阵运算结果相等。这让我们能构建线性方程组 $C_g c = y$ 来解出未知的 $c$，而且这个方程组同样可以用 FFT 快速求解。

### 证明

基于卷积定义（滑动窗口）和循环矩阵（每一行是上一行的移位）的概念，我们可以推导出循环矩阵为何能被傅里叶变换对角化。

推导的核心逻辑在于：**循环移位（Cyclic Shift）在傅里叶变换下变成了相位旋转（Phase Rotation）。**

#### 第一步：定义基本组件

1.  **循环矩阵 $C$**：
    由向量 $c = [c_0, c_1, \dots, c_{N-1}]^T$ 生成。矩阵的每一行都是上一行向右移动一位。
    我们可以把循环矩阵 $C$ 分解为**移位矩阵 $S$** 的多项式组合。

2.  **移位矩阵 $S$**：
    这是一个特殊的循环矩阵，它只把向量向右移动 1 位：
    $$
    S = \begin{bmatrix} 
    0 & 0 & \dots & 1 \\ 
    1 & 0 & \dots & 0 \\ 
    0 & 1 & \dots & 0 \\ 
    \vdots & \vdots & \ddots & \vdots 
    \end{bmatrix}
    $$
    
    * $S x$ 的结果是把 $x$ 循环右移 1 位。
    * $S^k x$ 的结果是把 $x$ 循环右移 $k$ 位。

3.  **用 $S$ 表示 $C$**：
    因为 $C$ 的第 $k$ 行（或其作用）等价于把卷积核移动 $k$ 步，所以循环矩阵可以写成：
    $$C = c_0 I + c_1 S + c_2 S^2 + \dots + c_{N-1} S^{N-1} = \sum_{k=0}^{N-1} c_k S^k$$

#### 第二步：寻找“移位”的特征向量

我们需要找到一个向量 $v$，使得它**移动一位后，仅仅是乘以了一个常数**（即 $Sv = \lambda v$）。

我们考察**傅里叶基向量**（Fourier Basis Vector）。设 $\omega = e^{-i \frac{2\pi}{N}}$ 是单位复根。
第 $k$ 个傅里叶基向量 $f_k$ 定义为：
$$
f_k = \begin{bmatrix} 1 \\ \omega^k \\ \omega^{2k} \\ \vdots \\ \omega^{(N-1)k} \end{bmatrix}
$$

如果我们把这个向量**向右循环移动一位**（即乘以 $S$），最下面的元素 $\omega^{(N-1)k}$ 会跑到最上面：

$$
S f_k = \begin{bmatrix} \omega^{(N-1)k} \\ 1 \\ \omega^k \\ \vdots \\ \omega^{(N-2)k} \end{bmatrix}
$$

利用性质 $\omega^{Nk} = 1$ (转一圈回到原点)，我们可以提取公因子 $\omega^{-k}$：

$$
S f_k = \omega^{-k} \begin{bmatrix} \omega^{Nk} \\ \omega^k \\ \omega^{2k} \\ \vdots \end{bmatrix} = \omega^{-k} \begin{bmatrix} 1 \\ \omega^k \\ \omega^{2k} \\ \vdots \end{bmatrix} = \omega^{-k} f_k
$$

**结论：**
傅里叶基向量 $f_k$ 是移位矩阵 $S$ 的**特征向量**！对应的特征值是 $\lambda_S = \omega^{-k}$。

#### 第三步：推导循环矩阵 $C$ 的对角化

既然 $C$ 只是 $S$ 的多项式组合（$C = \sum c_j S^j$），而 $f_k$ 是 $S$ 的特征向量，那么 **$f_k$ 必然也是 $C$ 的特征向量**。

计算 $C$ 乘以 $f_k$：

$$
\begin{aligned}
C f_k &= (\sum_{j=0}^{N-1} c_j S^j) f_k \\
&= \sum_{j=0}^{N-1} c_j (S^j f_k) \\
&= \sum_{j=0}^{N-1} c_j (\lambda_S)^j f_k  \quad \text{(代入S的特征值)} \\
&= \left( \sum_{j=0}^{N-1} c_j \omega^{-kj} \right) f_k
\end{aligned}
$$

**观察括号里的部分**：
$\sum_{j=0}^{N-1} c_j \omega^{-kj}$ 正好是向量 $c$ 的**离散傅里叶变换 (DFT)** 的第 $k$ 个分量！我们记作 $\hat{c}_k$。

所以我们得到了特征值方程：
$$C f_k = \hat{c}_k f_k$$

这意味着：
1.  **特征向量**：循环矩阵 $C$ 的特征向量就是傅里叶基向量 $f_k$。
2.  **特征值**：对应的特征值就是卷积核 $c$ 的傅里叶变换系数 $\hat{c}_k$。

#### 第四步：写成矩阵形式 (卷积定理)

我们将所有特征向量 $f_k$ 组合成一个矩阵，这正是**逆 DFT 矩阵** $F^{-1}$。
将所有特征值 $\hat{c}_k$ 放在对角线上，组成对角矩阵 $\Lambda = \text{diag}(\hat{c})$。

根据线性代数的特征分解公式 ($A = V \Lambda V^{-1}$)，我们要找的 $V$ 就是 $F^{-1}$。
即：
$$C = F^{-1} \cdot \text{diag}(\hat{c}) \cdot F$$

#### 最终结论总结

1.  **卷积运算** $y = c * x$ 等价于矩阵运算 $y = C x$。
2.  将 $C$ 替换为分解形式：
    $$y = F^{-1} \text{diag}(\hat{c}) F x$$
3.  两边同时左乘 $F$（做傅里叶变换）：
    $$F y = \text{diag}(\hat{c}) (F x)$$
4.  这等价于频域公式：
    $$\hat{y} = \hat{c} \odot \hat{x}$$

这就是为什么循环矩阵可以用 FFT 快速计算的数学本质：**因为循环矩阵在傅里叶基底下是对角的。**