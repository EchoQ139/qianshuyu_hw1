#include "algebra.h"
#include <stdio.h>
#include <math.h>

Matrix create_matrix(int row, int col)
{
    Matrix m;
    m.rows = row;
    m.cols = col;
    return m;
}

Matrix add_matrix(Matrix a, Matrix b)
{
    if(a.cols!=b.cols||a.rows!=b.rows) 
    {
        printf("Error: Matrix a and b must have the same rows and cols.\n");
        return create_matrix(0,0);
    }  
    Matrix result=create_matrix(a.rows,a.cols);
    for(int i=0;i<a.rows;i++)
    {
        for(int j=0;j<a.cols;j++)
            result.data[i][j]=a.data[i][j]+b.data[i][j];
    }
    return result;
}

Matrix sub_matrix(Matrix a, Matrix b)
{
    if(a.cols!=b.cols||a.rows!=b.rows) 
    {
        printf("Error: Matrix a and b must have the same rows and cols.\n");
        return create_matrix(0,0);
    }  
    Matrix result=create_matrix(a.rows,a.cols);
    for(int i=0;i<a.rows;i++)
    {
        for(int j=0;j<a.cols;j++)
            result.data[i][j]=a.data[i][j]-b.data[i][j];
    }
    return result;
}

Matrix mul_matrix(Matrix a, Matrix b)
{
    if(a.cols!=b.rows)
    {
        printf("Error: The number of cols of matrix a must be equal to the number of rows of matrix b.\n");
        return create_matrix(0, 0);
    }
    Matrix result=create_matrix(a.rows,b.cols);
    for(int i=0;i<a.rows;i++)
    {
        for(int j=0;j<b.cols;j++)
        {
            result.data[i][j]=0;
            for(int k=0;k<a.cols;k++)
            {
                result.data[i][j]+=a.data[i][k]*b.data[k][j];
            }
        }
    }
    return result;
}

Matrix scale_matrix(Matrix a, double k)
{
    Matrix result=create_matrix(a.rows,a.cols);
    for(int i=0;i<a.rows;i++)
    {
        for(int j=0;j<a.cols;j++)
            result.data[i][j]=a.data[i][j]*k;
    }
    return result;
}

Matrix transpose_matrix(Matrix a)
{
    Matrix result=create_matrix(a.cols,a.rows);
    for(int i=0;i<a.rows;i++)
    {
        for(int j=0;j<a.cols;j++)
            result.data[j][i]=a.data[i][j];
    }
    return result;
}

double calculate(Matrix a)
{
    double det=0;
    if(a.rows==1)
    {
        det=a.data[0][0];
        return det;
    }
    if(a.rows==2)
    {
        det=a.data[0][0]*a.data[1][1]-a.data[0][1]*a.data[1][0];
        return det;
    }
    for(int k=0;k<a.rows;k++)
    {
        Matrix sub_matrix=create_matrix(a.rows-1,a.cols-1);
        int subi=0;
        for(int i=0;i<sub_matrix.rows;i++)
        {
            int subj=0;
            if(i==k)
                continue;
            for(int j=1;j<sub_matrix.cols;j++)
            {
                sub_matrix.data[subi][subj]=a.data[i][j];
                subj++;
            }
            subi++;
        }
        det+=a.data[k][0]*pow(-1,k)*calculate(sub_matrix);
    }
    return det;
}

double det_matrix(Matrix a)
{
    if(a.rows!=a.cols)
    {
        printf("Error: The matrix must be a square matrix.\n");
        return 0;
    }
    return calculate(a);
}

Matrix inv_matrix(Matrix a)
{
    if(a.rows!=a.cols)
    {
        printf("Error: The matrix must be a square matrix.\n");
        return create_matrix(0, 0);
    }
    if(det_matrix(a)==0)
    {
        printf("Error: The matrix is singular.\n");
        return create_matrix(0, 0);
    }
    int n=a.rows;
    Matrix inv=create_matrix(n,n);
    double det=calculate(a);
    if(n==1)
    {
        inv.data[0][0]=1.0/a.data[0][0];
        return inv;
    }
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            Matrix sub_matrix=create_matrix(n-1,n-1);
            int subi=0;
            for(int q=0;q<n;q++)
            {
                if(q==i)
                    continue;
                int subj=0;
                for(int p=0;p<n;p++)
                {
                    if(p==j)
                        continue;
                    sub_matrix.data[subi][subj]=a.data[q][p];
                    subj++;
                }
                subi++;
            }
            double det_submatrix=pow(-1,i+j)*calculate(sub_matrix);
            inv.data[j][i]=det_submatrix/det;
        }
    }
    return inv;
}

Matrix swap_rows(Matrix a,int i1,int i2)
{
    double temp;
    for(int j=0;j<a.cols;j++)
    {
        temp=a.data[i1][j];
        a.data[i1][j]=a.data[i2][j];
        a.data[i2][j]=temp;
    }
}

int rank_matrix(Matrix a)
{
    int m=a.rows,n=a.cols;
    int rank;
    Matrix upper_triangle=create_matrix(m,n);
    for(int i=0;i<m;i++)
    {
        for(int j=0;j<n;j++)
            upper_triangle.data[i][j]=a.data[i][j];
    }
    if(m<=n)
        rank=m;
    else
        rank=n;
    for(int k=0;k<m&&k<n;k++)
    {
        if(upper_triangle.data[k][k]==0)
        {
            int flag=0;
            for(int q=k;q<m;q++)
            {
                if(upper_triangle.data[q][k]!=0)
                {
                    flag=1;
                    upper_triangle=swap_rows(upper_triangle,k,q);
                    break;
                }
            }
            if(flag==0)
            {
                rank--;
                
                continue;
            }
        }
        for(int i=0;i<m;i++)
        {
            if(i==k)
                continue;
            for(int j=k+1;j<n;j++)
            {
                 upper_triangle.data[i][j]-=upper_triangle.data[k][j]*upper_triangle.data[i][k]/upper_triangle.data[k][k];
            }
        }
        for(int i=0;i<m;i++)
            upper_triangle.data[i][k]=0;
    }
    return rank;
}

double trace_matrix(Matrix a)
{
    if(a.rows!=a.cols)
    {
        printf("Error: The matrix must be a square matrix.\n");
        return 0;
    }
    double trace=0;
    for(int i=0;i<a.rows;i++)
    {
        trace+=a.data[i][i];
    }
    return trace;
}

void print_matrix(Matrix a)
{
    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.cols; j++)
        {
            // 按行打印，每个元素占8个字符的宽度，小数点后保留2位，左对齐
            printf("%-8.2f", a.data[i][j]);
        }
        printf("\n");
    }
}