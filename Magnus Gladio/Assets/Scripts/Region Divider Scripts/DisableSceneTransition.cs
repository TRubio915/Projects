using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class DisableSceneTransition : MonoBehaviour
{
    private void DisableTransition()
    {
        this.gameObject.GetComponent<Animator>().enabled = false;
        this.gameObject.GetComponent<SpriteRenderer>().color = new Color(0, 0, 0, 1);
        this.gameObject.GetComponent<SpriteRenderer>().enabled = false;      

    }
}
