{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 19,
   "id": "b506aace",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "<matplotlib.image.AxesImage at 0x7f89e536e850>"
      ]
     },
     "execution_count": 19,
     "metadata": {},
     "output_type": "execute_result"
    },
    {
     "data": {
      "image/png": "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\n",
      "text/plain": [
       "<Figure size 432x288 with 1 Axes>"
      ]
     },
     "metadata": {
      "needs_background": "light"
     },
     "output_type": "display_data"
    }
   ],
   "source": [
    "import numpy as np\n",
    "import matplotlib.pyplot as plt\n",
    "my_data = np.genfromtxt('mass_grid_2d.csv', delimiter=',')\n",
    "my_data = my_data[:, :-1]\n",
    "\n",
    "\n",
    "\n",
    "plt.imshow(np.log(my_data + 1))"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "836baf84",
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3 (ipykernel)",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.9.5"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}